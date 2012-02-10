#include "ControlLogger.h"

void ControlLogger::ini() {
  latexFileName = new char[500];
  snprintf(latexFileName,500,"%s","logFile");
  loggingActive = false;
  plotCounter = 0;
  PLOT_PointSize = 1.5;
  PLOT_FrameSizeX = 0.8;
  PLOT_FrameSizeY = 0.8;
}


void ControlLogger::desini() {
  delete[] latexFileName;
  plotCounter = 0;
}


ControlLogger::ControlLogger() {
  ini();
}


ControlLogger::~ControlLogger() {
  desini();
}


void ControlLogger::addSection(const char* caption) {
  if (!loggingActive) return;
  char* fName = new char[1000];
  snprintf(fName,1000,"%s_content.tex",latexFileName);
  FILE* file = fopen(fName,"a");

  fprintf(file,"\\section*{%s}\n", caption);
  fclose(file);
  delete[] fName;
}


void ControlLogger::clearPage() {
  if (!loggingActive) return;
  char* fName = new char[1000];
  snprintf(fName,1000,"%s_content.tex",latexFileName);
  FILE* file = fopen(fName,"a");

  fprintf(file,"\\clearpage\n");
  fclose(file);
  delete[] fName;
}


void ControlLogger::newPage() {
  if (!loggingActive) return;
  char* fName = new char[1000];
  snprintf(fName,1000,"%s_content.tex",latexFileName);
  FILE* file = fopen(fName,"a");

  fprintf(file,"\\newpage\n");
  fclose(file);
  delete[] fName;
}



void ControlLogger::generateLatex() {
  if (!loggingActive) return;

  char* fName = new char[1000];
  snprintf(fName,1000,"%s.tex",latexFileName);
  FILE* file = fopen(fName,"w");

  fprintf(file,"\\documentclass[german,english,a4paper]{article}\n");
  fprintf(file,"\\usepackage[ngerman, english]{babel}\n");	
  fprintf(file,"\\usepackage[latin1]{inputenc}\n");
  fprintf(file,"\\usepackage[T1]{fontenc}\n");
  fprintf(file,"\\usepackage{ae}\n");	
  fprintf(file,"\\usepackage{amssymb}\n");
  fprintf(file,"\\usepackage{graphicx}\n");
  fprintf(file,"\\usepackage[figuresright]{rotating}\n");
  fprintf(file,"\\usepackage{subfigure}\n");
  fprintf(file,"\\usepackage{bbm}\n");
//  fprintf(file,"\\textwidth 155mm\n");
//  fprintf(file,"\\textheight 248mm\n");
//  fprintf(file,"\\topmargin -2.0cm \n");
  fprintf(file,"\\addtolength{\\textheight}{5.0cm}\n");
  fprintf(file,"\\addtolength{\\textwidth}{7.0cm}\n");
  fprintf(file,"\\addtolength{\\voffset}{-3.0cm}\n");
  fprintf(file,"\\addtolength{\\hoffset}{-3.0cm}\n");
  
  fprintf(file,"\\begin{document}\n");
  fprintf(file,"\\selectlanguage{english}\n");
  fprintf(file,"\\input{");
  fprintf(file,latexFileName);
  fprintf(file,"_content.tex}\n");
  fprintf(file,"\\end{document}\n");

  fclose(file);
  
  char* command = new char[1000];
  snprintf(command,1000,"latex %s",fName);
  
  system(command);
  delete[] fName;
  delete[] command;
}


void ControlLogger::setLatexFileName(char* fname) {
  if (!loggingActive) return;
  snprintf(latexFileName,1000,"%s",fname);
  char* command = new char[1000];
  snprintf(command,1000,"rm %s_content.tex",latexFileName);
  system(command);
  delete[] command;  
}


void ControlLogger::setLogging(bool on) {
  loggingActive = on;
}

  
void ControlLogger::addPlot(const char* title1, const char* title2, const char* caption, const char* xlabel, const char* ylabel, double* x, double* y1, double* err1, double* y2, double* err2, int N, const char* fitCommand) {
  if (!loggingActive) return;
  bool dob = true;
  if ((y2==NULL) || (err2==NULL) || (title2==NULL)) {
    y2 = y1;
    err2 = err1;
    title2 = title1;
    dob = false;
  }

  char* fName = new char[1000];
  snprintf(fName,1000,"%s%d%s","data/dataFile",plotCounter,".dat");
  FILE* file = fopen(fName,"w");
  
  int I;
  for (I=0; I<N; I++) {
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f %1.15f\n",x[I],y1[I],err1[I],y2[I],err2[I]);
  }
  fclose(file);
  
  snprintf(fName,1000,"pics/%s_plot_%d.gnu",latexFileName,plotCounter);
  file = fopen(fName,"w");
  fprintf(file,"set xlabel '%s'\n",xlabel);
  fprintf(file,"set ylabel '%s'\n",ylabel);
  fprintf(file,"set size %f,%f\n",PLOT_FrameSizeX,PLOT_FrameSizeY);    
  fprintf(file,"set pointsize %f\n",PLOT_PointSize);    
  fprintf(file,"plot 'data/dataFile%d.dat' using 1:2:3 with errorbars pt 6 title '%s'\n",plotCounter,title1);
  if (dob) {
    fprintf(file,"replot 'data/dataFile%d.dat' using 1:4:5 with errorbars pt 7 title '%s'\n",plotCounter,title2);
  }
  fprintf(file,fitCommand);
  fprintf(file,"\n");
  fclose(file);
  
  char* command = new char[1000];
  snprintf(command,1000,"./makePlot -bw pics/%s_plot_%d\n",latexFileName,plotCounter);
  system(command);
  delete[] command;

  fName = new char[1000];
  snprintf(fName,1000,"%s_content.tex",latexFileName);
  file = fopen(fName,"a");

  fprintf(file,"\\begin{figure}[h]\n");
  fprintf(file,"\\centering\n");
  fprintf(file,"\\includegraphics[angle=0,width=0.44\\textwidth]{pics/%s_plot_%d}\n",latexFileName,plotCounter);
  fprintf(file,"\\vspace{-5mm}\n");
  fprintf(file,"\\caption{%s}\n",caption);
  fprintf(file,"\\end{figure}\n");
  fprintf(file,"\\vspace{-1mm}\n\n");

  fclose(file);
  delete[] fName;
  plotCounter++;
}


void ControlLogger::addGnuplotTableLinePlot(char* title1, char* title2, char* caption, char* xlabel, char* ylabel, double* x, double* y1, double* y2, int N, char* fitCommand1, char* fitCommand2) {
  if (!loggingActive) return;

  char* fName = new char[1000];
  snprintf(fName,1000,"%s%d%s","data/dataFile",plotCounter,".dat");
  FILE* file = fopen(fName,"w");
  
  int I;
  for (I=0; I<N; I++) {
    fprintf(file,"%1.15f %1.15f %1.15f\n",x[I],y1[I],y2[I]);
  }
  fclose(file);
  
  snprintf(fName,1000,"pics/%s_plot_%d.gnu",latexFileName,plotCounter);
  file = fopen(fName,"w");
  fprintf(file,"set xlabel '%s'\n",xlabel);
  fprintf(file,"set ylabel '%s'\n",ylabel);
  fprintf(file,"set size %f,%f\n",PLOT_FrameSizeX,PLOT_FrameSizeY);    
  fprintf(file,"set pointsize %f\n",PLOT_PointSize);    
  fprintf(file,"set terminal postscript\n");    
  fprintf(file,"set output 'pics/%s_plot_%d.eps'\n",latexFileName,plotCounter);    
  fprintf(file,"plot 'data/dataFile%d.dat' using 1:2 with lines lt 6 title '%s'\n",plotCounter,title1);
  fprintf(file,"%s\n",fitCommand1);    
  fprintf(file,"set output 'pics/%s_plot_%d.eps'\n",latexFileName,plotCounter);    
  fprintf(file,"replot\n");    
  fprintf(file,"set output 'pics/%s_plot_%d.eps'\n",latexFileName,plotCounter+1);    
  fprintf(file,"plot 'data/dataFile%d.dat' using 1:3 with lines lt 6 title '%s'\n",plotCounter,title2);
  fprintf(file,"%s\n",fitCommand2);    
  fprintf(file,"set output 'pics/%s_plot_%d.eps'\n",latexFileName,plotCounter+1);    
  fprintf(file,"replot\n");    
  
  fprintf(file,"\n");
  fclose(file);
  
//  system("cd pics");
  char* command = new char[1000];
  snprintf(command,1000,"gnuplot pics/%s_plot_%d.gnu",latexFileName,plotCounter);
  system(command);
  delete[] command;
//  system("cd ..");

  fName = new char[1000];
  snprintf(fName,1000,"%s_content.tex",latexFileName);
  file = fopen(fName,"a");

  fprintf(file,"\\begin{figure}[h]\n");
  fprintf(file,"\\centering\n");
  fprintf(file,"\\begin{tabular}{cc}\n");
  fprintf(file,"\\includegraphics[angle=270,width=0.44\\textwidth]{pics/%s_plot_%d}  & \\includegraphics[angle=270,width=0.44\\textwidth]{pics/%s_plot_%d}\n",latexFileName,plotCounter,latexFileName,plotCounter+1);
  fprintf(file,"\\end{tabular}\n");
  
  fprintf(file,"\\vspace{-1mm}\n");
  fprintf(file,"\\caption{%s}\n",caption);
  fprintf(file,"\\end{figure}\n");
  fprintf(file,"\\vspace{-3mm}\n\n");

  fclose(file);
  delete[] fName;
  plotCounter+=2;
}


void ControlLogger::makePlotWithErrors(char* &fName, int N, double* x, double* y, double* err, char* xlabel, char* ylabel, char* title,char* fitCommand) {
  if (!loggingActive) return;

  fName = new char[1000];
  snprintf(fName,1000,"%s%d%s","data/dataFile",plotCounter,".dat");
  FILE* file = fopen(fName,"w");
  
  int I;
  for (I=0; I<N; I++) {
    fprintf(file,"%1.15f %1.15f %1.15f\n",x[I],y[I],err[I]);
  }
  fclose(file);
  
  snprintf(fName,1000,"pics/%s_plot_%d.gnu",latexFileName,plotCounter);
  file = fopen(fName,"w");
  fprintf(file,"set xlabel '%s'\n",xlabel);
  fprintf(file,"set ylabel '%s'\n",ylabel);
  fprintf(file,"set size %f,%f\n",PLOT_FrameSizeX,PLOT_FrameSizeY);    
  fprintf(file,"set pointsize %f\n",PLOT_PointSize);    
  fprintf(file,"plot 'data/dataFile%d.dat' using 1:2:3 with errorbars pt 6 title '%s'\n",plotCounter,title);
  fprintf(file,fitCommand);
  fprintf(file,"\n");
  fclose(file);
  
  char* command = new char[1000];
  snprintf(command,1000,"./makePlot -bw pics/%s_plot_%d",latexFileName,plotCounter);
  system(command);
  delete[] command;
  snprintf(fName,1000,"pics/%s_plot_%d",latexFileName,plotCounter);
  plotCounter++;
}
  
  
void ControlLogger::makePlotWithErrorsAndAdditionalLinePlots(char* &fName, int N, double* x, double* y1, double* err1, double* y2, double* err2, double* l1, double* l2, char* xlabel, char* ylabel, char* title1, char* title2, char* fitCommand) {
  if (!loggingActive) return;

  fName = new char[1000];
  snprintf(fName,1000,"%s%d%s","data/dataFile",plotCounter,".dat");
  FILE* file = fopen(fName,"w");
  
  int I;
  for (I=0; I<N; I++) {
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f\n",x[I],y1[I],err1[I],y2[I],err2[I],l1[I],l2[I]);
  }
  fclose(file);
  
  snprintf(fName,1000,"pics/%s_plot_%d.gnu",latexFileName,plotCounter);
  file = fopen(fName,"w");
  fprintf(file,"set xlabel '%s'\n",xlabel);
  fprintf(file,"set ylabel '%s'\n",ylabel);
  fprintf(file,"set size %f,%f\n",PLOT_FrameSizeX,PLOT_FrameSizeY);    
  fprintf(file,"set pointsize %f\n",PLOT_PointSize);    
  fprintf(file,"plot 'data/dataFile%d.dat' using 1:2:3 with errorbars pt 6 title '%s'\n",plotCounter,title1);
  fprintf(file,"replot 'data/dataFile%d.dat' using 1:4:5 with errorbars pt 7 title '%s'\n",plotCounter,title2);
  fprintf(file,"replot 'data/dataFile%d.dat' using 1:6 with lines lt 1 notitle\n",plotCounter);
  fprintf(file,"replot 'data/dataFile%d.dat' using 1:7 with lines lt 1 notitle\n",plotCounter);
  fprintf(file,fitCommand);
  fprintf(file,"\n");
  fclose(file);
  
  char* command = new char[1000];
  snprintf(command,1000,"./makePlot -bw pics/%s_plot_%d",latexFileName,plotCounter);
  system(command);
  delete[] command;
  snprintf(fName,1000,"pics/%s_plot_%d",latexFileName,plotCounter);
  plotCounter++;
}


void ControlLogger::addPlotTOTEX(char* eps1, char* caption, int angle1, double w1) {
  if (!loggingActive) return;

  char* fName = new char[1000];
  snprintf(fName,1000,"%s_content.tex",latexFileName);
  FILE* file = fopen(fName,"a");

  fprintf(file,"\\begin{figure}[h]\n");
  fprintf(file,"\\centering\n");
  fprintf(file,"\\includegraphics[angle=%d,width=%f\\textwidth]{%s} \n",angle1, w1, eps1);
  
  fprintf(file,"\\vspace{-3mm}\n");
  fprintf(file,"\\caption{%s}\n",caption);
  fprintf(file,"\\end{figure}\n");
  fprintf(file,"\\vspace{-1mm}\n\n");

  fclose(file);
  delete[] fName;
}

 
void ControlLogger::addDoublePlotTOTEX(char* eps1, char* eps2, char* caption, int angle1, double w1, int angle2, double w2) {
  if (!loggingActive) return;

  char* fName = new char[1000];
  snprintf(fName,1000,"%s_content.tex",latexFileName);
  FILE* file = fopen(fName,"a");

  fprintf(file,"\\begin{figure}[h]\n");
  fprintf(file,"\\centering\n");
  fprintf(file,"\\begin{tabular}{cc}\n");
  fprintf(file,"\\includegraphics[angle=%d,width=%f\\textwidth]{%s}  &  \\includegraphics[angle=%d,width=%f\\textwidth]{%s}\n",angle1, w1, eps1,angle2, w2,eps2);
  fprintf(file,"\\end{tabular}\n");
  
  fprintf(file,"\\vspace{-3mm}\n");
  fprintf(file,"\\caption{%s}\n",caption);
  fprintf(file,"\\end{figure}\n");
  fprintf(file,"\\vspace{-1mm}\n\n");

  fclose(file);
  delete[] fName;
}


void ControlLogger::addTablePlot(char* title1, char* title2, char* caption, char* xlabel1, char* xlabel2, char* ylabel1, char* ylabel2, double* x1, double* x2, double* y1, double* y2, double* yErr1, double* yErr2, int N1, int N2, char* fitCommand1, char* fitCommand2) {
  if (!loggingActive) return;
  char* eps1;
  char* eps2;

  makePlotWithErrors(eps1, N1, x1, y1, yErr1, xlabel1, ylabel1, title1, fitCommand1);
  makePlotWithErrors(eps2, N2, x2, y2, yErr2, xlabel2, ylabel2, title2, fitCommand2);
  
  addDoublePlotTOTEX(eps1, eps2, caption, 0, 0.44, 0, 0.44);

  delete[] eps1;
  delete[] eps2;
}


void ControlLogger::addPlotWithAdditionalLinePlots(char* title1, char* title2, char* caption, char* xlabel, char* ylabel, double* x, double* y1, double* err1, double* y2, double* err2, double* l1, double* l2, int N, char* fitCommand1, char* fitCommand2) {
  if (!loggingActive) return;
  char* eps1;

  makePlotWithErrorsAndAdditionalLinePlots(eps1, N, x, y1, err1, y2, err2, l1, l2, xlabel, ylabel, title1, title2, fitCommand1);
  addPlotTOTEX(eps1, caption, 0, 0.44);

  delete[] eps1;
}


void ControlLogger::addDirectText(char* text) {
  char* fName = new char[1000];
  snprintf(fName,1000,"%s_content.tex",latexFileName);
  FILE* file = fopen(fName,"a");

  fprintf(file,"%s\n",text);

  fclose(file);
  delete[] fName;
}


void ControlLogger::setPLOT_PointSize(double size) {
  PLOT_PointSize = size;
}


void ControlLogger::setPLOT_FrameSize(double x, double y) {
  PLOT_FrameSizeX = x;
  PLOT_FrameSizeY = y;
}
