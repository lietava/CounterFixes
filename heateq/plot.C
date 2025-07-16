{
	float x[8] = {1,10,50,100,200,300,500,1000};
	float y[8] = {39.1,33.1,26.4,24.0,22.8,22.64,22.6,22.6};
	TGraph *g = new TGraph(8,x,y);
	g->SetMarkerStyle(21);
	g->Draw("AP");
}
