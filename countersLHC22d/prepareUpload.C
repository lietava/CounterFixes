void prepareUpload() {

  FILE *fptr = fopen("command.txt", "w");

  std::vector<int> run = {527349, 527963, 528537, 528543};
  std::vector<long> startRoman = {1665784953893, 1666631792235, 1667374389991, 1667377507182};
  std::vector<long> startBK = {1665785484000, 1666632188000, 1667374928000, 1667378045000};
  std::vector<long> startTRG = {1665785553000, 1666632392000, 1667374990000, 1667378107000};
  std::vector<long> startGRP = {1665785484420, 1666632188522, 1667374928226, 1667378045403};

  std::vector<long> endBK = {1665786301000, 1666636722000, 1667377523000, 1667379469000};
  std::vector<long> endTRG = {1665786302000, 1666636722000, 1667377523000, 1667379469000};
  std::vector<long> endGRP = {1665786301934, 1666636722332, 1667377523129, 1667379469250};
  
  for (int i = 0; i < 4; ++i) {
    long startBKnew = startBK[i] - 1000 * 60 * 3;
    long startTRGnew = startTRG[i] - 1000 * 60 * 3;
    long startGRPnew = startGRP[i] - 1000 * 60 * 3;
    long endBKnew = endBK[i] + 1000 * 60 * 3;
    long endTRGnew = endTRG[i] + 1000 * 60 * 3;
    long endGRPnew = endGRP[i] + 1000 * 60 * 3;
    printf("\nrun %d\n", run[i]);
    printf("startBK    = %ld, startTRG    = %ld, startGRP    = %ld, endBK    = %ld, endRCT    = %ld, endGRP    = %ld\n", startBK[i], startTRG[i], startGRP[i], endBK[i], endTRG[i], endGRP[i]);
    printf("startBKnew = %ld, startTRGnew = %ld, startGRPnew = %ld, endBKnew = %ld, endTRGnew = %ld, endGRPnew = %ld\n", startBKnew, startTRGnew, startGRPnew, endBKnew, endTRGnew, endGRPnew);

    long start, end;
    start = std::min({startBKnew, startTRGnew, startGRPnew});
    end = std::max({endBKnew, endTRGnew, endGRPnew});

    printf("start = %ld end = %ld\n", start, end);

    fprintf(fptr, "o2-ccdb-upload -f %d.root --starttimestamp %ld --endtimestamp %ld  -k \"CTPRunScalers\" --path CTP/Calib/Scalers --host alice-ccdb.cern.ch -m \"JIRA=O2-3684;runNumber=%d\"\n", run[i], start, end, run[i]);
  }
  fclose(fptr);

}
