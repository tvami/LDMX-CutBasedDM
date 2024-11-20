import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python3 %prog sample.txt 3")
(opt,args) = parser.parse_args()

sampleInFile = sys.argv[1]
binIn = sys.argv[2]

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

cutFlowLin = True
cutFlowNorm = True

SamplesArray = []

with open(sampleInFile, "r") as a_file:
  for line in a_file:
    stripped_line = line.strip()
    SamplesArray.append(stripped_line)

fileInArray = []
for sample in SamplesArray:
  fileInArray.append(ROOT.TFile.Open(sample))
  
#dirs = []

#PostCutHistoBin = 13
PostCutHistoBin = int(binIn)
f = fileInArray[0]

for i in range(0, fileInArray[0].GetListOfKeys().GetEntries()):
  dirname = f.GetListOfKeys().At(i).GetName()
  curr_dir = f.GetDirectory(dirname)
#  print("dirname: "+dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      # Match the plot of interest
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      # if not ("CutFlow" in keyname):  continue
      keyname2 = keyname
      curr_dir2 = f.GetDirectory(dirname+"/"+keyname)
      if True :
          keyname = curr_dir.GetListOfKeys().At(i).GetName()
#          print(dirname+"/"+keyname)
          curr_dir2 = fileInArray[0].GetDirectory(dirname+"/"+keyname)
          newname = dirname+"/"+keyname
          obj = fileInArray[0].Get(newname)
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
          
          tex5 = ROOT.TLatex()
#          tex5 = ROOT.TLatex(0.07,0.04,"Code version: "+codeVersion);

          overFlowText = ROOT.TLatex(0.836,0.15,"OF");
#          overFlowText.SetTextAngle(-30)
          overFlowText.SetNDC();
          overFlowText.SetTextFont(52);
          overFlowText.SetTextSize(0.01);
          overFlowText.SetLineWidth(2);
          
          if obj.InheritsFrom("TObject"):
              if not os.path.exists(os.path.dirname("Compare"+sampleInFile[:-4]+"_Bin" + binIn + "/a.png")):
                print("Create dir")
                os.makedirs(os.path.dirname("Compare"+sampleInFile[:-4]+"_Bin" + binIn + "/"))
              if (obj.GetEntries() == 0 ) : continue

              # Now do the 2D histos
              if (obj.ClassName() == "TH2F"):
                canvasString = 'canvas'+str(i)
                canvas = ROOT.TCanvas(canvasString, canvasString, 800,800)
                if not ("TrigEff" in keyname or ("CutFlow" in keyname and cutFlowLin)) : canvas.SetLogy()
                #TADA
                legendStartX = 0.0
                legendEntryForZero = SamplesArray[0][0:SamplesArray[0].find(".root")]
                if ((len(legendEntryForZero)) < 25) :
                  legendStartX = 0.6
                if ((len(legendEntryForZero)) >= 25) :
                  legendStartX = 0.5
                if ((len(legendEntryForZero)) >= 40) :
                  legendStartX = 0.3
                if ((len(legendEntryForZero)) >= 55) :
                  legendStartX = 0.2
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
                
                histoArray = []
                # fake index to satisfy ROOT memory allocation
                i = 0
                for fileIn in fileInArray:
#                  print(keyname)
                  if ("N1_" in keyname) :
                    PostCutHisto = fileIn.Get(newname).ProjectionY(keyname2 + "_ProjY"+str(i),PostCutHistoBin,PostCutHistoBin)
                  elif ("Rev_" in keyname) :
                    PostCutHisto = fileIn.Get(newname).ProjectionY(keyname2 + "_ProjY"+str(i),PostCutHistoBin,PostCutHistoBin)
                  elif  ("TrigEff" in keyname):
                    PostCutHisto = fileIn.Get(newname).ProfileY(keyname2 + "_ProfY"+str(i))
                    for indexBin in range(1,PostCutHisto.GetNbinsX()):
                        if (PostCutHisto.GetBinError(indexBin) > 0.2 or PostCutHisto.GetBinError(indexBin) == 0 ) : PostCutHisto.SetBinContent(indexBin, -1.)
                  elif ("CutFlow" in keyname):
                    PostCutHisto = fileIn.Get(newname).ProjectionX(keyname2 + "_ProjX"+str(i))
                    PostCutHisto.LabelsOption("v")
                    # PostCutHisto.GetXaxis().SetBinLabel(2,"Non-fiducial")

                  else :
                    PostCutHisto = fileIn.Get(newname).ProjectionY(keyname2 + "_ProjY"+str(i),PostCutHistoBin,PostCutHistoBin)
                  
                  
                  if (PostCutHisto.Integral()> 0 and not "TrigEff" in keyname and not "CutFlow" in keyname) :
                    PostCutHisto.Scale(1/PostCutHisto.Integral(1,PostCutHisto.GetNbinsX()+1))
                  if (PostCutHisto.GetBinContent(1)>0 and "CutFlow" in keyname and cutFlowNorm) : 
                    PostCutHisto.Scale(1/PostCutHisto.GetBinContent(1))
                  if (PostCutHisto.GetBinContent(1)>0 and "CutFlow" in keyname) : 
                    print(f"| Cutflow | |")
                    print(f"| --- | --- | ")
                    for binIndex in range(0,PostCutHisto.GetNbinsX()+1):
                      bin_content = PostCutHisto.GetBinContent(binIndex)
                      bin_label = PostCutHisto.GetXaxis().GetBinLabel(binIndex)
                      if bin_content == 0.0: continue
                      print(f"| {bin_label} | {bin_content} |")

                  i += 1
                  if (PostCutHisto) : histoArray.append(PostCutHisto)
                for index in range(0, len(histoArray)):
                  histoArray[index].SetStats(0)
                  histoArray[index].SetMarkerStyle(20)
                  histoArray[index].GetXaxis().SetRange(1,histoArray[index].GetNbinsX()+1)
                  if ("PT" in keyname) :
                    histoArray[index].GetXaxis().SetRangeUser(0.0,1000.0)

                  legendEntry = SamplesArray[index][0:SamplesArray[index].find(".root")]
                  legend.AddEntry(histoArray[index],legendEntry,"LP")
                  indexNew = -1
                  if (index>-1):
                    indexNew = index+2
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
                  histoArray[index].SetLineColor(indexNew)
                  histoArray[index].SetMarkerColor(indexNew)
#                  histoArray[index].SetFillColor(indexNew)
                  histoArray[index].SetTitle("")
                  max = 0.0
#                  print(keyname2)
                  for index2 in range(0, len(histoArray)):
                    if not (histoArray[index2]) : continue
                    max = numpy.maximum(max,histoArray[index2].GetMaximum())
                  histoArray[0].GetYaxis().SetTitle("Normalized events / bin")
#                  histoArray[0].GetYaxis().SetRangeUser(0,max*1.05)
                  histoArray[0].GetYaxis().SetRangeUser(0.000001,max*5)
                  if ("TrigEff" in keyname) : 
                    histoArray[0].GetYaxis().SetRangeUser(0.0,1.5)
                    histoArray[0].GetYaxis().SetTitle("Efficiency")
                  if ("CutFlow" in keyname and cutFlowLin) : 
                    histoArray[0].GetYaxis().SetRangeUser(0.0,1.10)
                  binWidth = histoArray[index].GetXaxis().GetBinWidth(1)
                  firstBinLocXCent = histoArray[index].GetXaxis().GetBinCenter(1)
                  overFlowLocXCent = histoArray[index].GetXaxis().GetBinCenter(histoArray[index].GetNbinsX()+ 1)
                  overFlowLocX = (overFlowLocXCent - (binWidth)/2)
                  if ("CutFlow" in keyname) : 
                    overFlowLocX = 99999
                  firstBinLoxX = abs(firstBinLocXCent - (binWidth)/2)
                  
                  overFlowLine = ROOT.TLine();
#                  overFlowLine.SetNDC();
                  overFlowLine.SetLineWidth(2);
                  overFlowLine.SetLineStyle(ROOT.kDashed);
                  if "NReadoutHits" in keyname:
                    histoArray[index].Rebin(2)
                  histoArray[index].Draw("SAMEPE")
                  if ("CutFlow" in keyname and not cutFlowLin):
                    histoArray[index].Draw("SAMEHISTOTEXT45")
                  
#                  if (any(substring in keyname2 for substring in ["mass", "pt"])):
#                    histoArray[index].Rebin(2)

                legend.Draw("SAME")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                overFlowText.Draw("SAME")
                overFlowLine.DrawLine(overFlowLocX,histoArray[0].GetMinimum(),overFlowLocX,histoArray[0].GetMaximum())
                canvas.SaveAs("Compare"+sampleInFile[:-4]+"_Bin" + binIn + "/"+keyname2+".png")
                canvas.SaveAs("Compare"+sampleInFile[:-4]+"_Bin" + binIn + "/"+keyname2+".C")

#                 ROOT.gStyle.SetPadRightMargin(.09)
#                 ROOT.gStyle.SetPadTopMargin(0.1);
#                 ROOT.gStyle.SetPadBottomMargin(0.20);
#                 ROOT.gStyle.SetPadLeftMargin(0.15);
#                 canvasSCutFlowtring = 'canvasCutFlow'+str(i)
#                 canvasCutFlow = ROOT.TCanvas(canvasSCutFlowtring, canvasSCutFlowtring, 800,800)
#                 canvasCutFlow.SetLogy()
#                 cutflow.LabelsOption("v")
#                 cutflow.SetStats(0)
#                 cutflow.Draw("HISTOTEXT45")
#                 cutflow.GetYaxis().SetTitle("Events / bin")
#                 cutflow.GetYaxis().SetRangeUser(1,cutflow.GetMaximum()*100)
# #                if ("Rev" in keyname) :
# ##                  cutflow.GetXaxis().SetBinLabel(5,"E_{cell,max} < 300")
# #                  cutflow.GetXaxis().SetBinLabel(12,"Non-fiducial")
# #                else:
# ##                  cutflow.GetXaxis().SetBinLabel(10,"E_{cell,max} < 300")
# #                  cutflow.GetXaxis().SetBinLabel(3,"Non-fiducial")
#                 tex5.Draw("SAME")
#                 tex2.Draw("SAME")
#                 tex3.Draw("SAME")
#                 tex4.Draw("SAME")

#                 canvasCutFlow.SaveAs("Compare"+sampleInFile[:-5]+"/CutFlow"+str(i)+".png")

#                 canvasCutFlowNstring = 'canvasCutFlowN'+str(i)
#                 canvasCutFlowN = ROOT.TCanvas(canvasCutFlowNstring, canvasCutFlowNstring, 800,800)
#                 cutflow.SetStats(0)
#                 if (cutflow.GetBinContent(1)>0) : cutflow.Scale(1/cutflow.GetBinContent(1))
#                 cutflow.Draw("HISTOTEXT45")
#                 cutflow.GetYaxis().SetTitle("Efficiency")
# #                if ("Rev" in keyname) :
# ##                  cutflow.GetXaxis().SetBinLabel(5,"E_{cell,max} < 300")
# #                  cutflow.GetXaxis().SetBinLabel(12,"Non-fiducial")
# #                else:
# ##                  cutflow.GetXaxis().SetBinLabel(10,"E_{cell,max} < 300")
# #                  cutflow.GetXaxis().SetBinLabel(3,"Non-fiducial")
#                 cutflow.LabelsOption("v")
#                 tex2.Draw("SAME")
#                 tex3.Draw("SAME")
#                 tex4.Draw("SAME")
#                 tex5.Draw("SAME")
#                 canvasCutFlowN.SaveAs("Compare"+sampleInFile[:-5]+"/CutFlowNorm"+str(i)+".png")
