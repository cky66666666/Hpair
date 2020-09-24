#include "iostream"
#include "TH1D.h"
#include "THStack.h"
#include "TFile.h"
#include "TStyle.h"

void draw()
{
    TFile f("hist.root");
    THStack *a = (THStack*) f.Get("InvMass");
    a->Draw("nostack plc");
    gPad->BuildLegend();
}