#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TTree.h"

struct Nucleon {
  double x;
  double y;
  double z;
  bool participant;
};

// Build one nucleus by placing A nucleons independently and uniformly inside a
// hard sphere of radius R.
std::vector<Nucleon> SampleNucleus(int atomicNumber,
                                   double radius,
                                   TRandom3 &rng)
{
  std::vector<Nucleon> nucleons;
  nucleons.reserve(atomicNumber);

  while ((int)nucleons.size() < atomicNumber) {
    const double r = radius * std::cbrt(rng.Uniform());
    const double costheta = rng.Uniform(-1.0, 1.0);
    const double sintheta = std::sqrt(1.0 - costheta * costheta);
    const double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());

    Nucleon candidate;
    candidate.x = r * sintheta * std::cos(phi);
    candidate.y = r * sintheta * std::sin(phi);
    candidate.z = r * costheta;
    candidate.participant = false;

    nucleons.push_back(candidate);
  }

  return nucleons;
}

// Move the full nucleus in x by dx. In the collision setup this puts the two
// nuclear centers at -b/2 and +b/2 in the transverse plane.
void ShiftNucleus(std::vector<Nucleon> &nucleons, double dx)
{
  for (auto &n : nucleons) {
    n.x += dx;
    n.participant = false;
  }
}

// Draw a negative-binomial random number with mean mu and width parameter k.
// This uses the Gamma-Poisson representation:
//   lambda ~ Gamma(k, mu / k), then n ~ Poisson(lambda).
int SampleNegativeBinomial(double mean, double k, std::mt19937_64 &rng)
{
  if (mean <= 0.0 || k <= 0.0) {
    return 0;
  }

  std::gamma_distribution<double> gamma(k, mean / k);
  const double lambda = gamma(rng);
  std::poisson_distribution<int> poisson(lambda);
  return poisson(rng);
}

void glauber_mc(int nEvents = 10000,
                int atomicNumber = 208,
                double nucleusRadius = 6.62,
                double sigmaNNmb = 70.0,
                double bMax = 20.0,
                double meanParticlesPerParticipant = 20.0,
                double nbdK = 1.0,
                bool useNegativeBinomial = true,
                const char *outFileName = "glauber_output.root",
                unsigned int seed = 1)
{
  TRandom3 rng(seed);
  std::mt19937_64 nbdRng(seed + 12345);

  // 1 mb = 0.1 fm^2. Black-disk condition: d_T < sqrt(sigma_NN / pi).
  const double sigmaNNfm2 = 0.1 * sigmaNNmb;
  const double collisionDistance2 = sigmaNNfm2 / TMath::Pi();

  TFile outFile(outFileName, "RECREATE");

  TH1D hB("hB", "Impact parameter;b (fm);events", 200, 0.0, bMax);
  TH1D hNpart("hNpart", "Participants;N_{part};events", 2 * atomicNumber + 1, -0.5,
              2.0 * atomicNumber + 0.5);
  TH1D hNcoll("hNcoll", "Binary collisions;N_{coll};events", 1200, -0.5, 1199.5);
  TH1D hNparticles("hNparticles", "Produced particles;N_{particles};events",
                   12000, -0.5, 11999.5);
  TH2D hNpartNcoll("hNpartNcoll", "N_{coll} vs N_{part};N_{part};N_{coll}",
                   2 * atomicNumber + 1, -0.5, 2.0 * atomicNumber + 0.5,
                   1200, -0.5, 1199.5);
  TH2D hNpartNparticles("hNpartNparticles",
                        "Produced particles vs N_{part};N_{part};N_{particles}",
                        2 * atomicNumber + 1, -0.5, 2.0 * atomicNumber + 0.5,
                        12000, -0.5, 11999.5);
  TH2D hBNparticles("hBNparticles",
                    "Produced particles vs impact parameter;b (fm);N_{particles}",
                    200, 0.0, bMax, 12000, -0.5, 11999.5);

  TTree eventTree("events", "Glauber Monte Carlo event observables");
  int event = 0;
  double b = 0.0;
  int nPart = 0;
  int nColl = 0;
  int nParticles = 0;
  int multiplicityModel = useNegativeBinomial ? 1 : 0;
  eventTree.Branch("event", &event, "event/I");
  eventTree.Branch("b", &b, "b/D");
  eventTree.Branch("npart", &nPart, "npart/I");
  eventTree.Branch("ncoll", &nColl, "ncoll/I");
  eventTree.Branch("nparticles", &nParticles, "nparticles/I");
  eventTree.Branch("multiplicityModel", &multiplicityModel, "multiplicityModel/I");

  for (event = 0; event < nEvents; ++event) {
    // Sample b with probability proportional to b db, i.e. uniform in area.
    b = bMax * std::sqrt(rng.Uniform());

    auto projectile = SampleNucleus(atomicNumber, nucleusRadius, rng);
    auto target = SampleNucleus(atomicNumber, nucleusRadius, rng);

    ShiftNucleus(projectile, -0.5 * b);
    ShiftNucleus(target, 0.5 * b);

    // Count binary nucleon-nucleon collisions in the transverse plane.
    // If two nucleons collide, both become participants.
    nColl = 0;
    for (auto &p : projectile) {
      for (auto &t : target) {
        const double dx = p.x - t.x;
        const double dy = p.y - t.y;
        if (dx * dx + dy * dy < collisionDistance2) {
          p.participant = true;
          t.participant = true;
          ++nColl;
        }
      }
    }

    // Count all nucleons that had at least one binary collision.
    nPart = 0;
    for (const auto &p : projectile) {
      nPart += p.participant ? 1 : 0;
    }
    for (const auto &t : target) {
      nPart += t.participant ? 1 : 0;
    }

    // Particle-production model:
    //   useNegativeBinomial = true: fluctuate each participant with NBD.
    //   useNegativeBinomial = false: fixed particles per participant.
    if (useNegativeBinomial) {
      nParticles = 0;
      for (int i = 0; i < nPart; ++i) {
        nParticles += SampleNegativeBinomial(meanParticlesPerParticipant, nbdK, nbdRng);
      }
    } else {
      nParticles = (int)std::round(meanParticlesPerParticipant * nPart);
    }

    hB.Fill(b);
    hNpart.Fill(nPart);
    hNcoll.Fill(nColl);
    hNparticles.Fill(nParticles);
    hNpartNcoll.Fill(nPart, nColl);
    hNpartNparticles.Fill(nPart, nParticles);
    hBNparticles.Fill(b, nParticles);
    eventTree.Fill();
  }

  hB.Write();
  hNpart.Write();
  hNcoll.Write();
  hNparticles.Write();
  hNpartNcoll.Write();
  hNpartNparticles.Write();
  hBNparticles.Write();
  eventTree.Write();
  outFile.Close();

  std::cout << "Wrote " << nEvents << " Pb+Pb-style Glauber events to "
            << outFileName << std::endl;
  std::cout << "Parameters: A=" << atomicNumber
            << ", R=" << nucleusRadius << " fm"
            << ", sigmaNN=" << sigmaNNmb << " mb"
            << ", bMax=" << bMax << " fm"
            << ", meanParticlesPerParticipant=" << meanParticlesPerParticipant
            << ", nbdK=" << nbdK
            << ", useNegativeBinomial=" << useNegativeBinomial << std::endl;
}
