
#include <stdio.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <string> 
#include <iostream>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <math.h>

#include "LatexAndPlottingSystem.h"


using namespace std;
 
struct Constraint {
	string obsName;
	string tagName;
	string constraintRel;
	double constraintValue;
};
 
struct Column {
	string obsName;
	string tagName;
	string description;
};
 

struct XMLplot {
	string xmlFileName;
	string plotName;
	string plotTitle;
	string plotCaption;
	bool xLogScale;
	bool yLogScale;
	bool yErrorBars;
	bool xErrorBars;
	string title;
	string xLabel;
	string yLabel;
	
	double xSize, ySize;
	double pointSize;
	double xRangeMin, xRangeMax, yRangeMin, yRangeMax;
	int pointType, lineType;
	int lineWidth;
	
	vector<Column> columns;
	vector<string> usingStrings;
	vector<string> usingTitles;	
	vector<string> directPlots;			
	vector<vector<Constraint> > runConstraints;			
};

struct SummaryElement {
	string obsName;
	string tagName;
	string description;
	double value;
	double error;
};
 

bool readPlotDescription(xmlNode* root, XMLplot& xmlPlot);
SummaryElement* readSummaryElement(xmlNode* root, const string& obsName, const string& tagName);
xmlNode* findNode(xmlNode* root, const string& name );
void collectSummaryFiles(vector<string>& result, XMLplot& xmlPlot);
void printSubTree(xmlNode* root, int depth);
bool getContent(xmlNode* node, const string& token, string& result);
void getSummaryFileNames(vector<string>& summaryFileVector);
int str2int (const string &str, int defaultValue);
double str2double(const string& str, double defaultValue);
void extractDataFromSummaryFiles(const vector<string>& fileNames, vector< vector<double> >& data, XMLplot& xmlPlot);
void getFileNameList(char* searchDir, char* searchString, char** &fileNames, int& fileCount);
void deleteFileNameList(char** &fileNames, int& fileCount);
void printSummaryElement(SummaryElement& se);
void printXMLplot(XMLplot& xmlPlot);
void printArrayOfChar(const char* title, char** array, int count);
void printStringVector(char* title, vector<string>& stringVector);
void printDataArray(const char* title, const vector<vector<double> >& data);
void createFolders(XMLplot& xmlPlot, string& baseDirName);
void getFileNameFromPath(const string& path, string &fileName);
bool getContent(xmlNode* root, string& result);


int main(int argc, char **argv) {
	
	if (argc != 2) {
		fprintf(stderr, "ERROR: Missing xml file.\n");
		fprintf(stderr, "Usage: makeXMLplot <xml_file.xml>\n");
		fprintf(stderr, "subfolder XMLplots will be generated if not already present.\n");
		return(1);
	}
	system("clear");
	/*
	* this initialize the library and check potential ABI mismatches
	* between the version it was compiled for and the actual shared
	* library used.
	*/
	LIBXML_TEST_VERSION

	xmlDocPtr doc; /* the resulting document tree */	
	doc = xmlReadFile(argv[1], NULL, 0);
	if (doc == NULL) {
		fprintf(stderr, "Failed to parse %s\n", argv[1]);		
		exit(0);
	}
	xmlNode *root_element = NULL;
	root_element = xmlDocGetRootElement(doc);
	
	bool res = false;
	XMLplot xmlPlot;
		xmlPlot.xmlFileName = argv[1];
				
	if (readPlotDescription(root_element, xmlPlot)) {		
		printXMLplot(xmlPlot);
		
		xmlFreeDoc(doc);
		xmlCleanupParser();
		xmlMemoryDump();
		
		vector<string> summaryFileNames;	
		collectSummaryFiles(summaryFileNames, xmlPlot);				
		printStringVector("Summary-Files", summaryFileNames);
	
		vector<vector<double> > data;
		extractDataFromSummaryFiles(summaryFileNames, data, xmlPlot);
//		printDataArray("Data-Array", data);
		
		string baseDirName;
		createFolders(xmlPlot, baseDirName);		

		LatexAndPlottingSystem lapSystem(baseDirName.c_str(), xmlPlot.plotName.c_str());
		string plotScriptPrefix("plot");
		LAPsystemPlot* plot = lapSystem.createNewPlot(plotScriptPrefix.c_str());
		
			int dataLineCount = data.size();
			int dataColumnCount = xmlPlot.columns.size()*2;
	
			double** plotData = new double*[dataLineCount];
			for (int i = 0; i < dataLineCount; i++) {
				plotData[i] = new double[dataColumnCount];
				for (int j = 0; j < dataColumnCount; j++) {				
					plotData[i][j] = data.at(i).at(j);
				}
			}		
			
			plot->setPlotData(dataLineCount, dataColumnCount, plotData);
					
			for (int i = 0; i < dataLineCount; i++) {
				delete[] plotData[i];
				plotData[i] = NULL;
			}
			delete[] plotData;
			plotData = NULL;
	
			plot->setXLabel(xmlPlot.xLabel.c_str());		
			plot->setYLabel(xmlPlot.yLabel.c_str());				 				
			plot->setCaption(xmlPlot.plotCaption.c_str());
			plot->setPlotTitle(xmlPlot.plotTitle.c_str());
			plot->setXLogScale(xmlPlot.xLogScale);
			plot->setYLogScale(xmlPlot.yLogScale);
			plot->setXRange(xmlPlot.xRangeMin, xmlPlot.xRangeMax);
			plot->setYRange(xmlPlot.yRangeMin, xmlPlot.yRangeMax);
			plot->setXErrorBars(xmlPlot.xErrorBars);
			plot->setYErrorBars(xmlPlot.yErrorBars);  
			plot->setSize(xmlPlot.xSize, xmlPlot.ySize);
			for(size_t i = 0; i < xmlPlot.usingStrings.size(); i++) {
				plot->setTitle(xmlPlot.usingTitles.at(i).c_str());
				plot->plotData(xmlPlot.usingStrings.at(i).c_str());
			}
			if(xmlPlot.directPlots.size() != 0) {
				for (size_t i = 0; i < xmlPlot.directPlots.size(); i++) {					
					plot->plotDirect(xmlPlot.directPlots.at(i).c_str());
				}				
			}
		lapSystem.addPlot(plot);
			
			
		res = true;
	}
	else {
		cout << "ERROR while reading plot description! Check XML file name" << endl;
		xmlFreeDoc(doc);
		xmlCleanupParser();
		xmlMemoryDump();
	}
	
	return(res);
}

bool readPlotDescription(xmlNode* root, XMLplot& xmlPlot) {
	string str;
	xmlNode* curr_node = NULL;
	curr_node = root->children;
	bool okay = true;
	
	xmlNode* node = NULL;
	
	node = findNode(curr_node, "plotName");	
	if (!getContent(root, "plotName", xmlPlot.plotName) )
		okay = false;
	if (xmlPlot.plotName.size() == 0) {
		cout << "ERROR in xml file. plotName is empty!" << endl;
		okay = false;
	}
	
	getContent(root, "plotTitle", xmlPlot.plotTitle);
	getContent(root, "plotCaption", xmlPlot.plotCaption);
		
	xmlNode* plotInfo = findNode(curr_node, "plotInfo");
		getContent(plotInfo, "xLogScale", str);
		xmlPlot.xLogScale = str2int(str, 0);
		
		getContent(plotInfo, "yLogScale", str);
		xmlPlot.yLogScale = str2int(str, 0);
				
		getContent(plotInfo, "xLabel", str);
		xmlPlot.xLabel = str;
		
		getContent(plotInfo, "yLabel", str);
		xmlPlot.yLabel = str;
		
		getContent(plotInfo, "yErrorBars", str);
		xmlPlot.yErrorBars = str2int(str, 0);
		
		getContent(plotInfo, "xErrorBars", str);
		xmlPlot.xErrorBars = str2int(str, 0);
		
		getContent(plotInfo, "xSize", str);
		xmlPlot.xSize = str2double(str, 1.0);
		
		getContent(plotInfo, "ySize", str);
		xmlPlot.ySize = str2double(str, 1.0);
		
		getContent(plotInfo, "pointSize", str);
		xmlPlot.pointSize = str2int(str, 1);
		
		getContent(plotInfo, "xRangeMin", str);
		xmlPlot.xRangeMin = str2double(str, -INFINITY);
		
		getContent(plotInfo, "xRangeMax", str);
		xmlPlot.xRangeMax = str2double(str, +INFINITY);
		
		getContent(plotInfo, "yRangeMin", str);
		xmlPlot.yRangeMin = str2double(str, -INFINITY);
		
		getContent(plotInfo, "yRangeMax", str);
		xmlPlot.yRangeMax = str2double(str, INFINITY);
		
		getContent(plotInfo, "pointType", str);
		xmlPlot.pointType = str2int(str, 1);
		
		getContent(plotInfo, "lineType", str);
		xmlPlot.lineType = str2int(str, 1);
		
		getContent(plotInfo, "lineWidth", str);
		xmlPlot.lineWidth= str2int(str, 1);
	
	node = findNode(root, "plotData");
	while (node != NULL) {
		getContent(node, "usingString", str);
		if (str.size() != 0) {			
			xmlPlot.usingStrings.push_back(str);
			getContent(node, "title", str);
			xmlPlot.usingTitles.push_back(str);
		}
		node = findNode(node->next, "plotData");
	}	
		
		
	node = findNode(root, "directPlot");
	while (node != NULL) {
		getContent(node, str);
		if (str.size() != 0) {
			xmlPlot.directPlots.push_back(str);
		}
		node = findNode(node->next, "directPlot");
	}
		
	
	Column col;
	node = findNode(root, "column");
	if (node == NULL)
			cout << "No column found!" << endl;
	while (node != NULL) {
		getContent(node, "obsName", col.obsName);
		getContent(node, "tagName", col.tagName);

		xmlPlot.columns.push_back(col);
		node = findNode(node->next, "column");
	}
		
	xmlNode* runSelector = findNode(root, "runSelector");
	xmlNode* cons_node = NULL;
	vector<Constraint> consVector;
	Constraint cons;
	while (runSelector != NULL) {		
		cons_node = findNode(runSelector, "constraint");
		consVector.clear();				
		while (cons_node != NULL) {
			getContent(cons_node, "obsName", cons.obsName);
			getContent(cons_node, "tagName", cons.tagName);
			getContent(cons_node, "constraintRel", cons.constraintRel);
			getContent(cons_node, "constraintValue", str);
			cons.constraintValue = str2double(str, 0.0);
			consVector.push_back(cons);							
			cons_node = findNode(cons_node->next, "constraint");
		}	
		xmlPlot.runConstraints.push_back(consVector);
		runSelector = findNode(runSelector->next, "runSelector");
	}
	
	/*
	Constraint cons;
	xmlNode* runSelector = findNode(curr_node, "runSelector");
	xmlNode* cons_node = findNode(runSelector, "constraint");
	while (cons_node != NULL) {
		getContent(cons_node, "obsName", cons.obsName);
		getContent(cons_node, "tagName", cons.tagName);
		getContent(cons_node, "constraintRel", cons.constraintRel);
		getContent(cons_node, "constraintValue", str);
		cons.constraintValue = str2double(str, 0.0);
			
		xmlPlot.constraints.push_back(cons);
		cons_node = findNode(cons_node->next, "constraint");
	}
	*/
	return okay;
}

SummaryElement* readSummaryElement(xmlNode* root, const string& obsName, const string& tagName) {
	string str;
	SummaryElement* se = NULL;
	
	if (root != NULL && obsName.size() != 0 && tagName.size() != 0) {
		xmlNode* node = findNode(root->children, obsName);
		if (node != NULL) {
			node = findNode(node->children, tagName);
			if (node != NULL) {
				se = new SummaryElement();
				se->obsName = obsName;
				se->tagName = tagName;
				
				str = "";
				getContent(node, "value", str);
if (se->obsName.compare("Lambda0") == 0)
	cout << "------------> " << str << endl;				
				se->value = str2double(str, 0.0);
				
				str = "";
				getContent(node, "error", str);
				se->error = str2double(str, 0.0);
				
				getContent(node, "description", se->description);
			}
		}
	}
	return se; 
}

bool getContent(xmlNode* node, const string& token, string& result) {	
	string str;
	result = "";
	bool okay = false;
	
	if (node != NULL && token.size() != 0) {
		xmlNode* curr_node = NULL;
		bool found = false;
		for (curr_node = node->children; curr_node; curr_node = curr_node->next) {
			if (curr_node->type == XML_ELEMENT_NODE) {
				str = (char*)curr_node->name;
				if (str == token) {
					char* tmp = (char*)xmlNodeGetContent(curr_node->children);
					if (tmp != NULL && strlen(tmp) != 0) {
						result = tmp;
						okay = true;
					}
					xmlFree(tmp);
					found = true;
					break;
				}		
			}	
		}	
	}
	return okay;
}

bool getContent(xmlNode* root, string& result) {	
	string str;
	result = "";
	bool okay = false;
	
	if (root != NULL) {
		xmlNode* curr_node = root;
		if (curr_node->type == XML_ELEMENT_NODE) {
			char* tmp = (char*)xmlNodeGetContent(curr_node->children);
			if (tmp != NULL && strlen(tmp) != 0) {
				result = tmp;
				okay = true;
			}
			xmlFree(tmp);
		}			
	}
	return okay;
}


xmlNode* findNode(xmlNode* root, const string& name ) {
	xmlNode* res = NULL;		
	xmlNode *cur_node = NULL;
	char* node_name;
	bool found = false;
	
	for (cur_node = root; cur_node && !found; cur_node = cur_node->next) {
		if (cur_node->type == XML_ELEMENT_NODE) {
			node_name = (char*)cur_node->name;
			if (strcmp(node_name, name.c_str()) == 0) {
				res = cur_node;
				found = true;
			}
		}
		if (res == NULL)
 				res = findNode(cur_node->children, name);
	}
	return res;
}

void collectSummaryFiles(vector<string>& result, XMLplot& xmlPlot) {
	
	vector<string> summaryFileNames;
	getSummaryFileNames(summaryFileNames); 			
	SummaryElement* se = NULL;
	xmlDocPtr summaryDoc; /* the resulting document tree */
	xmlNode* summary_root = NULL;
	
	bool constraintFulfilled = false;
	
	for (size_t i = 0; i < summaryFileNames.size(); i++) {					
		summaryDoc = xmlReadFile(summaryFileNames.at(i).c_str(), NULL, 0);
		if (summaryDoc == NULL) {
			fprintf(stderr, "Failed to parse %s\n", summaryFileNames.at(i).c_str());
			exit(0);
		}		
		summary_root = xmlDocGetRootElement(summaryDoc);
		
		for(size_t runs = 0; runs < xmlPlot.runConstraints.size(); runs++) {
			constraintFulfilled = true;
			vector<Constraint> constraints(xmlPlot.runConstraints.at(runs));
			for (size_t j = 0; j < constraints.size(); j++) {				
				se = readSummaryElement(summary_root, constraints.at(j).obsName, constraints.at(j).tagName);						
				if (se != NULL) {
					//printSummaryElement(*se);
					//check if constraints are fulfilled				
					if (constraints.at(j).constraintRel.compare("=") == 0) {
						if (isnan(constraints.at(j).constraintValue) || isinf(constraints.at(j).constraintValue)) {
							if ( !isnan(se->value) ) {
								constraintFulfilled = false;
							}
						}
						else if ( (fabs(constraints.at(j).constraintValue - se->value) > se->error) ) {
							constraintFulfilled = false;			
						}					
					}				
					else if (constraints.at(j).constraintRel.compare("<") == 0 || constraints.at(j).constraintRel.compare("le") == 0) {
						if ( !((se->value - constraints.at(j).constraintValue) <= se->error) ) {
							constraintFulfilled = false;	
						}
					}
					else if (constraints.at(j).constraintRel.compare(">") == 0 || constraints.at(j).constraintRel.compare("ge") == 0 ) {
						if ( !((se->value - constraints.at(j).constraintValue) + se->error >= 0.0) ) {
							constraintFulfilled = false;
						}
					}					
				}			
				else {
					constraintFulfilled = false;
				}
								
				delete se;
				se = NULL;	
			}
			if (constraintFulfilled)
				result.push_back(summaryFileNames.at(i));		
		}		
				
			
		xmlFreeDoc(summaryDoc);
		xmlCleanupParser();
		xmlMemoryDump();
	}
	
	
}

void printSubTree(xmlNode* root, int depth) {
	xmlNode* curr_node = NULL;
	string str;
	
	string whiteSpace;
	for (int i = 0; i < depth; i++) {
		whiteSpace +="\t";
	}
	
	for (curr_node = root; curr_node; curr_node = curr_node->next) {
		if (curr_node->type == XML_ELEMENT_NODE) {
			str = (char*)curr_node->name;
			cout << whiteSpace << "=> " << str << endl;
		}
		printSubTree(curr_node->children, depth + 1);
	}			
}

int str2int (const string &str, int defaultValue) {	
	istringstream iss (str, istringstream::in);
	int n;
	if (!(iss >> n))
		n = defaultValue;
	return n;
}

double str2double(const string& str, double defaultValue) {
	double value = 0.0;
	if (str.size() != 0) {
		int res = sscanf (str.c_str(),"%lf", &value); 
		if (res == 0 || res == EOF) 
			value = defaultValue;
	
		/*	
		istringstream iss (str, istringstream::in);	
		if(!(iss >> value))
			value = defaultValue;
		*/
	}
	return value;
}

void getSummaryFileNames(vector<string>& summaryFileVector) {
	char* searchDir = new char[2048];
	char* searchString = new char[2048];
	char** fileNames = NULL;
	int fileCount = 0;
	
	sprintf(searchDir, "dataBase/data/results/pHMC/analysis");
	sprintf(searchString, "subFolderL");
	
	getFileNameList(searchDir, searchString, fileNames, fileCount);
	//printArrayOfChar("Available Sub-Folders", fileNames, fileCount);
	char* subSearchDir = new char[2048];
	char* subSearchString = new char[2048];
	char** summaryFiles = NULL;
	int summaryFileCount = 0;
	for (int i = 0; i < fileCount; i++ ) {			
		sprintf(subSearchDir, "%s", fileNames[i]);
		sprintf(subSearchString, "SummaryL");
	
		getFileNameList(subSearchDir, subSearchString, summaryFiles, summaryFileCount);
		for (int j = 0; j < summaryFileCount; j++) {
			summaryFileVector.push_back(summaryFiles[j]);
		}
		deleteFileNameList(summaryFiles, summaryFileCount);			
	}	
	delete[] subSearchDir;
	delete[] subSearchString;
	
	deleteFileNameList(fileNames, fileCount);
	delete[] searchDir;
	delete[] searchString;	
}

void extractDataFromSummaryFiles(const vector<string>& fileNames, vector< vector<double> >& data, XMLplot& xmlPlot) {
	SummaryElement* se = NULL;
	xmlDocPtr summaryDoc; /* the resulting document tree */
	xmlNode* summary_root = NULL;
	data.clear();
	vector<double> dataLine(xmlPlot.columns.size()*2, 0.0);
	int posDataColCount = 0;
	for (size_t i = 0; i < fileNames.size(); i++) {
		summaryDoc = xmlReadFile(fileNames.at(i).c_str(), NULL, 0);
		if (summaryDoc == NULL) {
			fprintf(stderr, "Failed to parse %s\n", fileNames.at(i).c_str());
			exit(0);
		}		
		summary_root = xmlDocGetRootElement(summaryDoc);
		
		posDataColCount = 0;
		for (size_t j = 0; j < xmlPlot.columns.size(); j++) {
			se = readSummaryElement(summary_root, xmlPlot.columns.at(j).obsName, xmlPlot.columns.at(j).tagName);						
			if (se != NULL) {
				dataLine.at(2*j) = se->value;
				dataLine.at(2*j+1) = se->error;
				posDataColCount++;				
			}
			else {
				dataLine.at(2*j) = 0.0;
				dataLine.at(2*j+1) = 0.0;
				string no_path_file_name("");
				getFileNameFromPath(fileNames.at(i), no_path_file_name);
				cout << "NOTICE: Could not read column " << xmlPlot.columns.at(j).obsName << "->" << xmlPlot.columns.at(j).tagName << endl;
				cout << "from file " << no_path_file_name << endl;
				cout << endl;
			}
			delete se;
			se = NULL;
		}
		if ((size_t)posDataColCount == xmlPlot.columns.size() && posDataColCount > 0)
			data.push_back(dataLine);
		//xmlPlot.columns.at(j).description = se->description;
		
		xmlFreeDoc(summaryDoc);
		xmlCleanupParser();
		xmlMemoryDump();
	}
	cout << "Non-Zero data lines: " << data.size() << endl;
}

void getFileNameList(char* searchDir, char* searchString, char** &fileNames, int& fileCount) {
	fileCount = 0;
	fileNames = NULL;

	for (int loop=0; loop<2; loop++) {
		DIR *dp = NULL;
		struct dirent *ep;
  
		if (searchDir == NULL) {
			dp = opendir("."); 
		} else {
			dp = opendir(searchDir);
		}

		if (dp == NULL) return;
  
		while ((ep = readdir (dp)) != NULL) {
			bool match = true;

			if (searchString != NULL) {
				if (strlen(ep->d_name) < strlen(searchString)) {
					match = false;
				} else {
					for (int I=0; I<((int)strlen(searchString)); I++) {
						if (searchString[I] != ep->d_name[I]) { 
							match = false;
							break;
						}
					}
				}
			}

			if (match) {
				if (loop==1) {
					fileNames[fileCount] = new char[strlen(ep->d_name)+2+strlen(searchDir)];
					snprintf(fileNames[fileCount], strlen(ep->d_name)+2+strlen(searchDir), "%s/%s", searchDir, ep->d_name);
				}
				fileCount++;
			}
		}
    
		if (loop==0) {
			fileNames = new char*[fileCount];
			fileCount = 0;
		}
    
		closedir(dp);
	}
}

void deleteFileNameList(char** &fileNames, int& fileCount) {
	for (int I=0; I<fileCount; I++) {
		delete[] fileNames[I];    
	}
	delete[] fileNames;
	fileNames = NULL;
	fileCount = 0;
}

void printSummaryElement(SummaryElement& se) {
	printf("Sunmmary element:\n");
	printf("%20s: %s\n", "Observablename", se.obsName.c_str());
	printf("%20s: %s\n", "Tagname", se.tagName.c_str());
	printf("%20s: %f (%f)\n", "Value", se.value, se.error);
	printf("%20s: %s\n", "Description", se.description.c_str());	
}

void printXMLplot(XMLplot& xmlPlot) {
	printf("XML Plot Description:\n");
	printf("%30s: %s\n", "Plotname", xmlPlot.plotName.c_str());
	printf("%30s: %s\n", "Plot Title", xmlPlot.plotTitle.c_str());
	printf("%30s: %s\n", "Plot Caption", xmlPlot.plotCaption.c_str());
	printf("%30s: %d\n", "Use log scale on x-Axis", xmlPlot.xLogScale);
	printf("%30s: %d\n", "Use log scale on y-Axis", xmlPlot.yLogScale);
	printf("%30s: %d\n", "Use error bars on x-Axis", xmlPlot.xErrorBars);
	printf("%30s: %d\n", "Use error bars on y-Axis", xmlPlot.yErrorBars);
	printf("%30s: %s\n", "x-Axis Label", xmlPlot.xLabel.c_str());
	printf("%30s: %s\n", "y-Axis Label", xmlPlot.yLabel.c_str());
	printf("%30s: (%f,%f)\n", "Size (x,y)", xmlPlot.xSize, xmlPlot.ySize);
	printf("%30s: %f\n", "Point Size", xmlPlot.pointSize);
	printf("%30s: (%f,%f)\n", "X-Range (min, max)", xmlPlot.xRangeMin, xmlPlot.xRangeMax);
	printf("%30s: (%f,%f)\n", "Y-Range (min, max)", xmlPlot.yRangeMin, xmlPlot.yRangeMax);
	printf("%30s: %d\n", "Line Type", xmlPlot.lineType);
	printf("%30s: %d\n", "Point Type", xmlPlot.pointType);
	printf("%30s: %d\n", "Line Width", xmlPlot.lineWidth);		
	
	printf("\n%30s:\n", "COLUMNS");
	for (size_t i = 0; i < xmlPlot.columns.size(); i++) {
		printf("%30s: %s\n", "Observable-Name", xmlPlot.columns.at(i).obsName.c_str());
		printf("%30s: %s\n\n", "Tag-Name", xmlPlot.columns.at(i).tagName.c_str());
	}
	
	printf("\n%30s:\n", "RUN-SELECTOR");
	for (size_t r = 0; r < xmlPlot.runConstraints.size(); r++) {
		vector<Constraint> constraints(xmlPlot.runConstraints.at(r));
		for (size_t i = 0; i < constraints.size(); i++) {
			printf("%30s: %s\n", "Observable-Name", constraints.at(i).obsName.c_str());
			printf("%30s: %s\n", "Tag-Name", constraints.at(i).tagName.c_str());
			printf("%30s: %s%f\n\n", "value", constraints.at(i).constraintRel.c_str(), constraints.at(i).constraintValue);
		}
		if ((r+1) < xmlPlot.runConstraints.size())
			printf("\n%30s:\n", "OR"); 
				
	}	
	
	printStringVector("Plot Column(s)", xmlPlot.usingStrings);
	printStringVector("Column Title", xmlPlot.usingTitles);
	printStringVector("Direct Plot(s)", xmlPlot.directPlots);
}

void printArrayOfChar(const char* title, char** array, int count) {
	printf("%s:\n", title);
	for (int i = 0; i < count; i++) {
		printf("%3d: %s\n", (i+1), array[i]);
	}
}

void printStringVector(char* title, vector<string>& stringVector) {
	printf("%3s:\n", title);
	for (size_t i = 0; i < stringVector.size(); i++) {
		printf("%3d.) %s\n", (int)(i+1), stringVector.at(i).c_str());
	}
}

void printDataArray(const char* title, const vector<vector<double> >& data) {
	printf("%s:\n", title);
	
	if(data.size() > 0 && data.at(0).size() > 0) {
		printf("%-10s\t", "Line no.");
		for (size_t i = 0; i < data.at(0).size(); i++) {
			printf("%-10s %d\t", "Column", (int)i);
		}
		printf("\n");
		
		for (size_t i = 0;  i < data.size(); i++) {
			printf("%5d%10s\t", (int)(i+1), "");
			for (size_t j = 0; j < data.at(i).size(); j++) {
				printf("%-1.10f\t", data.at(i).at(j));
			}
			printf("\n");
		}
	}
	
}

void createFolders(XMLplot& xmlPlot, string& baseDirName) {
	string fileNamePrefix;
	getFileNameFromPath(xmlPlot.xmlFileName, fileNamePrefix);
		
	baseDirName = "./XMLplots/";	
	baseDirName.append(fileNamePrefix);
	cout << "Base directory name: " << baseDirName << endl;	
	
	string sysCmd("mkdir ./XMLplots");
	system(sysCmd.c_str());
	
	sysCmd = "mkdir ";
	sysCmd.append(baseDirName);		
	system(sysCmd.c_str());
	
	sysCmd = "mkdir ";
	sysCmd.append(baseDirName);
	sysCmd.append("/plots");		
	system(sysCmd.c_str());
	
	sysCmd = "cp ./makePlot ";
	sysCmd.append(baseDirName);
	sysCmd.append("/plots");		
	system(sysCmd.c_str());
}

void getFileNameFromPath(const string& path, string &fileName) {
	size_t endIndex, startIndex;
	startIndex = path.find_last_of("/");
	if (startIndex == string::npos)
		startIndex = 0;
	else 
		startIndex++;
		
	endIndex = path.find_last_of(".");	
	if (endIndex == string::npos ) {
		cout << "Fatal error: could not extract file name prefix from path!" << endl;
		cout << "Path: " << path <<  " Start index: " << startIndex << " End index: " << endIndex << endl;
		cout << path.substr( startIndex, endIndex) << endl;
	}
	else
		endIndex -= startIndex;	
	
	fileName = path.substr( startIndex, endIndex);
}
