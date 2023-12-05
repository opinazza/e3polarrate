void inspect(const char* tel = "POLA-01", bool extracut=1){
  TChain *t = new TChain("stats");
  char input[200];
  t->AddFile(Form("out%s.root",tel));

  int nev = t->GetEntries();
  printf("nev = %d\n",nev);
  // pair mapping
  // pair_index    channel_indexes_in_pair
  // pair=0  - bottom L=7,R=2   - top L=5,R=0
  // pair=1  - bottom L=10,R=15 - top L=5,R=0
  // pair=2  - bottom L=11,R=6  - top L=5,R=0
  // pair=3  - bottom L=14,R=3  - top L=5,R=0
  // pair=4  - bottom L=7,R=2   - top L=8,R=13
  // pair=5  - bottom L=10,R=15 - top L=8,R=13
  // pair=6  - bottom L=11,R=6  - top L=8,R=13
  // pair=7  - bottom L=14,R=3  - top L=8,R=13
  // pair=8  - bottom L=7,R=2   - top L=9,R=4
  // pair=9  - bottom L=10,R=15 - top L=9,R=4
  // pair=10 - bottom L=11,R=6  - top L=9,R=4
  // pair=11 - bottom L=14,R=3  - top L=9,R=4
  // pair=12 - bottom L=7,R=2   - top L=12,R=1
  // pair=13 - bottom L=10,R=15 - top L=12,R=1
  // pair=14 - bottom L=11,R=6  - top L=12,R=1
  // pair=15 - bottom L=14,R=3  - top L=12,R=1
  int ch1[16] = {7, 10, 11, 14, 7, 10, 11, 14, 7, 10, 11, 14, 7, 10, 11, 14}; // index bottom left vs pair_index
  int ch2[16] = {2, 15, 6, 3, 2, 15, 6, 3, 2, 15, 6, 3, 2, 15, 6, 3};         // index bottom right vs pair_index
  int ch3[16] = {5, 5, 5, 5, 8, 8, 8, 8, 9, 9, 9, 9, 12, 12, 12, 12};         // index top left vs pair_index
  int ch4[16] = {0, 0, 0, 0, 13, 13, 13, 13, 4, 4, 4, 4, 1, 1, 1, 1};         // index top right vs pair_index

  int VertPairCh[16] = {0, 3, 0, 3, 2, 0, 2, 0, 1, 2, 1, 2, 3, 1, 3, 1}; // vertical pair_index corresponding to each channel

  // define variables to point tree leafs (each tree entry corresponds to 1 minutes)
  double ts; t->SetBranchAddress("time",&ts);                 // timestamp from 1 Jan 2007
  int status; t->SetBranchAddress("status",&status);          // status (0=good minute)
  float duration; t->SetBranchAddress("duration",&duration);  // sampling duration inside the minute (60 = all seconds were acquired) 
  float rateRaw; t->SetBranchAddress("rateRaw",&rateRaw);     // raw rate (=trigger rate)  (majority condition)
  float rate; t->SetBranchAddress("rate",&rate);                      // rate majority condition + 1 single track
  float rate4c; t->SetBranchAddress("rate4c",&rate4c);                // rate requiring all 4 SiPM fired + 1 single track
  float ratePair[16]; t->SetBranchAddress("ratePair",ratePair);       // rate for each pair of plates (majority condition + 1 single track)
  float ratePair4c[16]; t->SetBranchAddress("ratePair4c",ratePair4c); //  "    "   "    "   "   "     (4AND condition + 1 single track) 
  float pres; t->SetBranchAddress("pres",&pres);                      // pressure in mbar
  float lat; t->SetBranchAddress("lat",&lat);                         // latitude
  float avTot[16]; t->SetBranchAddress("avTot",avTot);                // average ToT in the current minute
  float eff[16]; t->SetBranchAddress("eff",eff);                      // pseudo-efficiency for each channel  4AND/majority
  float parExtra[3]; t->SetBranchAddress("parExtra",parExtra);        // extra parameters from slow control output
  float parRates[2]; t->SetBranchAddress("parRates",parRates);        // rates from slot control output
  
  t->GetEvent(0); // get first time in the tree
  double start = ts;
  t->GetEvent(nev-1); // get last time in the tree
  double stop = ts;
  printf("%lf to %lf\n",start,stop);
  int hours = (stop - start)/3600/24;

  // define histos
  
  TH2F *hRateVsLat = new TH2F("hRateVsLat",";latitude (#circ); rate (Hz)",900,0,90,100,25,45);
  TProfile *hRateVsLatProj = new TProfile("hRateVsLatProy",";latitude (#circ); rate (Hz)",900,0,90);
  TProfile *hRawTimeNC = new TProfile("hRawTimeNC",Form("%s; date; rate (Hz)",tel),hours,start,stop);
  TProfile *hRawTime = new TProfile("hRawTime",Form("%s; date; rate (Hz)",tel),hours,start,stop);
  TProfile *hTimeMajCorr = new TProfile("hTimeMajCorr",Form("%s; date; rate (Hz)",tel),hours,start,stop);
  //  TH2F *hTimeLat = new TH2F("hTimeLat",Form("%s; date; latitude",tel),hours,start,stop,180,0,90);
  TProfile *hTimeLat = new TProfile("hTimeLat",Form("%s; date; latitude",tel),hours,start,stop);
  TProfile *hTime4cCorr = new TProfile("hTime4cCorr",Form("%s; date; rate (Hz)",tel),hours,start,stop);
  hRawTimeNC->SetLineColor(6);
  hTimeMajCorr->SetLineColor(2);
  hTime4cCorr->SetLineColor(4); 
  
  TH2F *hRateRaw = new TH2F("hRateRaw",";rate slow; rate reco",500,0,50,500,0,50);
  TH1F *hRateRawDiff = new TH1F("hRateRawDiff",";rate slow - rate reco (s)",100,-10,10);

  TH1F *hcEff = new TH1F("hcEff","; eff;",100,0,1);
  
  TProfile *hEff[16];
  TProfile *hEffMaj[16];
  for(int ipair=0; ipair < 16; ipair++){
    hEff[ipair] = new TProfile(Form("hEff_%d",ipair),Form("Pair %d; eff_4C; rate4CNorm",ipair),100,0,1);
    hEffMaj[ipair] = new TProfile(Form("hEffMaj_%d",ipair),Form("Pair %d; eff_maj; rateNorm",ipair),100,0.9,1);
    hEff[ipair]->SetLineWidth(3);
    hEffMaj[ipair]->SetLineWidth(3);
  }

  // define parameters
  float barCoef = 2.2E-3;
  float presRef = 1000;
  
  float refVert[16], ref4vert, refVertRef[4];
  for(int i=0; i < nev; i++){
    t->GetEvent(i);
    if(pres < 10) continue; // discard event without a valid value for pressure

    float presCorr = TMath::Exp(barCoef*(pres-presRef));

    // fill rate correcte per pressure vs time for all events with a valid value for pressure
    hRawTimeNC->Fill(ts, rateRaw*presCorr);

    if(!status) continue; // remove bad quality events

    if(fabs(rateRaw - parRates[0]) > 2) continue; // require raw rate from reconstruction consistent with the one in the slow control (slow control value not perfectly in sync with reco) within 2 Hz

    // fill raw rate
    hRateRaw->Fill(parRates[0], rateRaw);

    // compare raw rate and the value stored in slow control output
    hRateRawDiff->Fill(parRates[0] - rateRaw);

    // fill raw rate corrected for pressure vs time after cuts applied
    hRawTime->Fill(ts, rateRaw*presCorr);


    // from now on let's try to correct rate also for channel efficiencies
    float corrRate = 0;
    float corrRate4c = 0;
    bool isWrong = false;

    ref4vert = 0;
    for(int ipair=0; ipair < 16; ipair++){ // loop over all pairs of plates and correct each single rate for efficiency -> sum all rates
      refVert[ipair] = 0;
      if(rate < 10 || ((rate -rate4c)/rate > 0.1 && extracut)) continue; // some requirement to accept the event: rate > 10 and 4AND vs Majority within 10% (if required)

      float ceff = eff[ch1[ipair]]*eff[ch2[ipair]]*eff[ch3[ipair]]*eff[ch4[ipair]]; // efficiency for 4AND
      if(ceff > 0.2){
	corrRate4c += ratePair4c[ipair]/ceff;
      } else {
	isWrong = true; // all pairs at least with efficiency > 20% 
      }

      // efficiency for Majority
      ceff = eff[ch1[ipair]]*eff[ch2[ipair]]*eff[ch3[ipair]];
      ceff += eff[ch1[ipair]]*eff[ch2[ipair]]*eff[ch4[ipair]];
      ceff += eff[ch1[ipair]]*eff[ch3[ipair]]*eff[ch4[ipair]];
      ceff += eff[ch2[ipair]]*eff[ch3[ipair]]*eff[ch4[ipair]];
      ceff -= 3*eff[ch1[ipair]]*eff[ch2[ipair]]*eff[ch3[ipair]]*eff[ch4[ipair]];
      if(ceff > 0.2){
	corrRate += ratePair[ipair]/ceff;
	if(ipair%5 == 0){
	  ref4vert += ratePair[ipair]/ceff;              // sum rate for all vertical pair (0, 5, 10, 15)
	  refVertRef[ipair/5] = ratePair[ipair]/ceff;    // keep rates for each vertical pair
	}
      } else {
	isWrong = true; // all pairs at least with efficiency > 20% 
      }
    }

    // this is for 
    for(int ipair=0; ipair < 16; ipair++){
      refVert[ipair] = ref4vert;
      int chP1 = ch1[ipair];
      int chP2 = ch3[ipair];
      if(VertPairCh[chP1] == VertPairCh[chP2]){ // remove rate of the correspondig pair to avoid bias (this is used as a refernce to check that efficiency correction are good)
	refVert[ipair] -= refVertRef[VertPairCh[chP1]];
	refVert[ipair] /= 3;	
      } else {                                  // remove rates of the corresponding pairs to avoid bias (this is used as a refernce to check that efficiency correction are good)
	refVert[ipair] -= refVertRef[VertPairCh[chP1]];
	refVert[ipair] -= refVertRef[VertPairCh[chP2]];
	refVert[ipair] /= 2;	
      }
    } // refVert[ipair] tells you what is the expected vertical rate for this specific ipair without using ipair-rate in the calculation 

    
    if(isWrong || corrRate < 15 || (rate/corrRate)<0.95) continue; // remove bad quality events in the comparison

    hcEff->Fill(rate/corrRate);
    
    for(int ipair=0; ipair < 16; ipair++){ // report not corrected rate vs efficiency once normalied to the expected rate from other pairs
      float ceff = eff[ch1[ipair]]*eff[ch2[ipair]]*eff[ch3[ipair]]*eff[ch4[ipair]];
      if(ceff > 0.5){
	hEff[ipair]->Fill(ceff, ratePair4c[ipair] / refVert[ipair]);
      }

      ceff = eff[ch1[ipair]]*eff[ch2[ipair]]*eff[ch3[ipair]];
      ceff += eff[ch1[ipair]]*eff[ch2[ipair]]*eff[ch4[ipair]];
      ceff += eff[ch1[ipair]]*eff[ch3[ipair]]*eff[ch4[ipair]];
      ceff += eff[ch2[ipair]]*eff[ch3[ipair]]*eff[ch4[ipair]];
      ceff -= 3*eff[ch1[ipair]]*eff[ch2[ipair]]*eff[ch3[ipair]]*eff[ch4[ipair]];
      if(ceff > 0.5){
	hEffMaj[ipair]->Fill(ceff, ratePair[ipair] / refVert[ipair]);
      }
    }
    
    hTimeMajCorr->Fill(ts, corrRate*presCorr);
    hTime4cCorr->Fill(ts, corrRate4c*presCorr);
    if(lat > 35){
      hTimeLat->Fill(ts, lat/1.4);

      if(ts > 497.5E6 && ts < 499.5E6){
	float rr = corrRate*presCorr;
	hRateVsLat->Fill(lat, rr);
	if(rr > 25 && rr < 45) hRateVsLatProj->Fill(lat, rr);
      }
    }
    
  }

  // draw and save outputs
  
  new TCanvas;
  hRateRaw->Draw("colz");
  new TCanvas;
  hRateRawDiff->Draw();

  TCanvas *ceff = new TCanvas;
  TF1 *ff = new TF1("ff","[0]*x",0,1);
  ceff->Divide(4,4);
  for(int ipair=0; ipair < 16; ipair++){
    ceff->cd(ipair+1);
    hEff[ipair]->Draw();
    hEff[ipair]->SetStats(0);
    hEff[ipair]->Fit(ff);
    hEff[ipair]->SetTitle(Form("Fraction = %.2f%c",ff->GetParameter(0)*100,'%'));
  }  

  TCanvas *ceffMaj = new TCanvas;
  ceffMaj->Divide(4,4);
  for(int ipair=0; ipair < 16; ipair++){
    ceffMaj->cd(ipair+1);
    hEffMaj[ipair]->Draw();
    hEffMaj[ipair]->SetStats(0);
    hEffMaj[ipair]->Fit(ff);
    hEffMaj[ipair]->SetTitle(Form("Fraction = %.2f%c",ff->GetParameter(0)*100,'%'));
  }

  new TCanvas;
  hRawTimeNC->Draw();
  hRawTimeNC->SetLineWidth(4);
  hRawTime->Draw("SAME");
  hRawTime->SetLineWidth(3);
  hRawTimeNC->SetStats(0);
  hRawTimeNC->GetXaxis()->SetTimeDisplay(1);
  hRawTimeNC->GetXaxis()->SetTimeFormat("%d/%m/%y%F2007-01-01 00:00:00s0");
  hTimeMajCorr->Draw("SAME");
  hTimeMajCorr->SetLineWidth(3);
  hTime4cCorr->Draw("SAME");

  new TCanvas;
  TH1D *hProRaw = hRawTime->ProjectionX();
  TH1D *hProMaj = hTimeMajCorr->ProjectionX();
  TH1D *hPro4c = hTime4cCorr->ProjectionX();
  hProMaj->Divide(hProMaj, hProRaw, 1, 1, "B");
  hPro4c->Divide(hPro4c, hProRaw, 1, 1, "B");
  hProMaj->Draw();
  hProMaj->GetXaxis()->SetTimeDisplay(1);
  hProMaj->GetXaxis()->SetTimeFormat("%d/%m/%y%F2007-01-01 00:00:00s0");
  hPro4c->Draw("same");
  hProMaj->SetLineColor(2);
  hPro4c->SetLineColor(4);

  new TCanvas;
  hcEff->Draw();
  TFile *fout = new TFile(Form("res%s.root",tel),"RECREATE");
  hTimeMajCorr->Write();
  hTimeLat->Write();
  hRateVsLat->Write();
  hRateVsLatProj->Write();
  fout->Close();
}
