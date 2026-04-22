float p2(float mu)
{
  float pexp = exp(-mu);
  float p2 = 1 - pexp - mu*pexp;
  return p2;
}
constexpr int n = 11;
void tail()
{
  //int run[] = {570128, 570143, 570158, 570159, 570160, 570162, 570163, 570165, 570166, 570191, 570205, 570243, 570279, 570292};
  //float mu[] = {0.041, 0.041,0. 032,   0.032,  0.032,  0.032,  0.032,  0.032,  0.031,  0.032,  0.032,  0.032,  0.032,  0.031 };
  // p2 = 1- p0 -p1
  float x[n], y[n];
  float mumin = 0.03;
  float mumax = 0.041;
  float step = (mumax -  mumin) / n;
  float mu = mumin;
  for(int i = 0; i < n; i++){
    x[i] = mu;
    y[i] = p2(mu);
    mu += step;
  }
  TGraph *g =  new TGraph(n,x,y);
  g->SetTitle("P > 1");
  g->GetXaxis()->SetTitle("#mu");
  g->SetMarkerSize(10);
  g->Draw("AP");
}
