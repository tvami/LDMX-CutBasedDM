import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python3 %prog sample.txt")
(opt,args) = parser.parse_args()

sampleInFile = sys.argv[1]

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetPalette(1)

SamplesArray = []

fileIn = ROOT.TFile.Open(sampleInFile)

for i in range(0, fileIn.GetListOfKeys().GetEntries()):
  dirname = fileIn.GetListOfKeys().At(i).GetName()
  curr_dir = fileIn.GetDirectory(dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      newname = dirname+"/"+keyname
      # if not ("Trig" in keyname) : continue
      if True :
          keyname2 = keyname
#          print(newname)
          ROOT.gStyle.SetPadRightMargin(.15)
          ROOT.gStyle.SetPadTopMargin(0.1);
          ROOT.gStyle.SetPadBottomMargin(0.14);
          ROOT.gStyle.SetPadLeftMargin(0.15);
          obj = fileIn.Get(newname)
          obj.SetMarkerStyle(20)
          
          tex2 = ROOT.TLatex(0.15,0.92,"LDMX");
          tex2.SetNDC();
          tex2.SetTextFont(61);
          tex2.SetTextSize(0.0675);
          tex2.SetLineWidth(2);

          
          tex3 = ROOT.TLatex(0.33,0.92,"Simulation"); # for square plots
          tex3.SetNDC();
          tex3.SetTextFont(52);
          tex3.SetTextSize(0.0485);
          tex3.SetLineWidth(2);
          
          tex4 = ROOT.TLatex(0.64,0.92,"1e14 EOT (8 GeV)")
          tex4.SetNDC();
          tex4.SetTextFont(52);
          tex4.SetTextSize(0.03);
          tex4.SetLineWidth(2);
          
          ratioLen = len(sampleInFile)+4
          startTex5 = (120-ratioLen)/120
          print(startTex5)
          tex5 = ROOT.TLatex(startTex5,0.015,"Sample: "+sampleInFile[:-5]);
          tex5.SetNDC()
          tex5.SetTextFont(52)
          tex5.SetTextSize(0.017)
          tex5.SetLineWidth(2)

          overFlowText = ROOT.TLatex(0.836,0.15,"OF");
#          overFlowText.SetTextAngle(-30)
          overFlowText.SetNDC();
          overFlowText.SetTextFont(52);
          overFlowText.SetTextSize(0.01);
          overFlowText.SetLineWidth(2);
          
          if obj.InheritsFrom("TObject"):
              if not os.path.exists(os.path.dirname("Compare"+sampleInFile[:-5]+"/a.png")):
                print("Create dir")
                os.makedirs(os.path.dirname("Compare"+sampleInFile[:-5]+"/"))
              if (obj.GetEntries() == 0 ) : continue

              # Now do the 2D histos
              if (obj.ClassName() == "TH2F"):
              
                canvasString = 'canvas'+str(i)
                canvas = ROOT.TCanvas(canvasString, canvasString, 800,800)
                canvas.SetLogy()
                legendStartX = 0.5
                legend =  ROOT.TLegend(legendStartX,.75,.85,.89,"","brNDC")
                legend.SetTextFont(42)
                legend.SetTextSize(0.017)
                legend.SetBorderSize(1);
                legend.SetBorderSize(0);
                legend.SetLineColor(1);
                legend.SetLineStyle(1);
                legend.SetLineWidth(1);
                legend.SetFillColor(0);
#                legend.SetFillStyle(1001);
                legend.SetFillStyle(0);
                
                myHist = fileIn.Get(newname)
                if not (myHist) : continue
                cutflow = myHist.ProjectionX()
                max = 0.0
                rangeMax = myHist.GetNbinsX() if "_Hcal_MaxSector" in keyname else myHist.GetNbinsY()
                for y in range(1,rangeMax+1) :
                  projY = myHist.ProjectionY(keyname2 + "_ProjY"+str(y),y,y)
                  
                  
                  
                  projY.SetMarkerStyle(20)
                  if projY.Integral()==0 : continue
                  legendEntry = "Cut#"+str(y-1) + ": " + str(projY.Integral(1,projY.GetNbinsX()+1)) + ", avg = " +  str(round(projY.GetMean(),2)) + "#pm" + str(round(projY.GetStdDev(),2))
                  projY.Scale(1/projY.Integral())
                  projY.SetStats(0)
                  projY.Draw("SAMEP")
                  projY.GetXaxis().SetRange(1,projY.GetNbinsX()+1)
                  if ("RecoilX" in keyname) :
                    projY.GetXaxis().SetRange(0,projY.GetNbinsX()+1)
                  
                  max = numpy.maximum(max,projY.GetMaximum()*1000)
                  projY.GetYaxis().SetRangeUser(0.0000001,max)
                  projY.GetYaxis().SetTitle("Normalized vents / bin")
                  legend.AddEntry(projY,legendEntry,"LP")
                  indexNew = -1
                  if (y>-1):
                    indexNew = y+1
#                  print(indexNew)
                  if (indexNew==10) :
                    indexNew = 40
                  elif (indexNew==11) :
                    indexNew = 46
                  elif (indexNew==12) :
                    indexNew = 41
                  elif (indexNew==13) :
                    indexNew = 30
                  elif (indexNew==14) :
                    indexNew = 42
#                  projY.SetFillColor(indexNew)
#                  projY.SetTitle("")
                  projY.SetLineColor(indexNew)
                  projY.SetMarkerColor(indexNew)
                    
                  binWidth = projY.GetXaxis().GetBinWidth(1)
                  firstBinLocXCent = projY.GetXaxis().GetBinCenter(1)
                  overFlowLocXCent = projY.GetXaxis().GetBinCenter(projY.GetNbinsX()+ 1)
                  overFlowLocX = (overFlowLocXCent - (binWidth)/2)
                  firstBinLoxX = abs(firstBinLocXCent - (binWidth)/2)
                  
                  overFlowLine = ROOT.TLine();
#                  overFlowLine.SetNDC();
                  overFlowLine.SetLineWidth(2);
                  overFlowLine.SetLineStyle(ROOT.kDashed);
#                  if (any(substring in keyname2 for substring in ["mass", "pt"])):

                legend.Draw("SAME")
                #tex2.Draw("SAME")
                #tex3.Draw("SAME")
                #tex4.Draw("SAME")
                tex5.Draw("SAME")
                overFlowText.Draw("SAME")
                overFlowLine.DrawLine(overFlowLocX,projY.GetMinimum(),overFlowLocX,projY.GetMaximum())
#                canvas.SetLogy()
#                print(keyname2)
                canvas.SaveAs("Compare"+sampleInFile[:-5]+"/"+keyname2+".png")
                
                ROOT.gStyle.SetPadRightMargin(.09)
                ROOT.gStyle.SetPadTopMargin(0.1);
                ROOT.gStyle.SetPadBottomMargin(0.20);
                ROOT.gStyle.SetPadLeftMargin(0.15);
                canvasSCutFlowtring = 'canvasCutFlow'+str(i)
                canvasCutFlow = ROOT.TCanvas(canvasSCutFlowtring, canvasSCutFlowtring, 800,800)
                canvasCutFlow.SetLogy()
                cutflow.LabelsOption("v")
                cutflow.SetStats(0)
                cutflow.Draw("HISTOTEXT45")
                cutflow.GetYaxis().SetTitle("Events / bin")
                cutflow.GetYaxis().SetRangeUser(1,cutflow.GetMaximum()*100)
#                if ("Rev" in keyname) :
##                  cutflow.GetXaxis().SetBinLabel(5,"E_{cell,max} < 300")
#                  cutflow.GetXaxis().SetBinLabel(12,"Non-fiducial")
#                else:
##                  cutflow.GetXaxis().SetBinLabel(10,"E_{cell,max} < 300")
                # cutflow.GetXaxis().SetBinLabel(2,"Non-fiducial")
                # cutflow.GetXaxis().SetBinLabel(3,"Triggerred")
                tex5.Draw("SAME")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")

                canvasCutFlow.SaveAs("Compare"+sampleInFile[:-5]+"/CutFlow"+str(i)+".png")

                canvasCutFlowNstring = 'canvasCutFlowN'+str(i)
                canvasCutFlowN = ROOT.TCanvas(canvasCutFlowNstring, canvasCutFlowNstring, 800,800)
                cutflow.SetStats(0)
                if (cutflow.GetBinContent(1)>0) : cutflow.Scale(1/cutflow.GetBinContent(1))
                cutflow.Draw("HISTOTEXT45")
                cutflow.GetYaxis().SetTitle("Efficiency")
#                if ("Rev" in keyname) :
##                  cutflow.GetXaxis().SetBinLabel(5,"E_{cell,max} < 300")
#                  cutflow.GetXaxis().SetBinLabel(12,"Non-fiducial")
#                else:
##                  cutflow.GetXaxis().SetBinLabel(10,"E_{cell,max} < 300")
#                  cutflow.GetXaxis().SetBinLabel(3,"Non-fiducial")
                cutflow.LabelsOption("v")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                canvasCutFlowN.SaveAs("Compare"+sampleInFile[:-5]+"/CutFlowNorm"+str(i)+".png")
