#include "gnuplotcall.h"

void replaceUnderscoreWithDot(string& str) {
    for (char& c : str) {
        if (c == '_') {
            c = '.';
        }
    }
}

string getContentAfterSlash(const string& inputString) {
    size_t slashPosition = inputString.find('/');
    if (slashPosition != string::npos) {
        return inputString.substr(slashPosition + 1);
    }
    else {
        return inputString;
    }
}

void ResultPresentation(vector<string> &sv_files1, string &outputPathname, vector<string> &tool_names, string &outputBasicMetricschart, vector< vector<float> > MeticsValues, vector< vector<int> > MeticsValues1){

	cout << endl <<"############# Stage 5: Generate comparison graph information: #############" << endl << endl;
	if(sv_files1.size()>1){
		string newInfo = "tool Recall Precision F1-score";
		string newInfo1 = "tool TP_user TP_benchmark FP FN";
		outputBasicMetricschartPath = outputPathname + outputBasicMetricschart;
		mkdir(outputBasicMetricschartPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		multipledataset(MeticsValues, sv_files1, tool_names, outputBasicMetricschartPath, newInfo);
		multipledataset(MeticsValues1, sv_files1, tool_names, outputBasicMetricschartPath, newInfo1);

		outputDiffRangeBasicMetricschart = outputBasicMetricschartPath + '/' + outputDiffRangeBasicMetricschart;
		mkdir(outputDiffRangeBasicMetricschart.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		//Draw multiple single graph patterns
//		ComparisonofMetricsindifferentrangesOp(outputDiffRangeBasicMetricschart, tool_names, sv_files1);
		//Draw a multi-graph mode
		vector<string> filenames;
		GenerateFileOp(outputDiffRangeBasicMetricschart, tool_names, filenames, sv_files1);
		GenerateMultiBarCharts(outputDiffRangeBasicMetricschart, filenames, newInfo);

		//sv information comparison of different sizes
		SVsizeratiofile = outputBasicMetricschartPath + '/' + SVsizeratiofile;
		mkdir(SVsizeratiofile.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		GenerateSVsizeRatioFileGraph(sv_files1, tool_names, SVsizeratiofile);
		//html
		Generatehtml(outputBasicMetricschart);
		cout << endl << "Details of the assessment results are saved in:" << outputPathname + htmlFilename<<endl;
	}else{
		outputBasicMetricschartPath = outputPathname + outputBasicMetricschart;
		mkdir(outputBasicMetricschartPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		Histogram_drawing(MeticsValues,  outputPathname, outputBasicMetricschartPath);
		Histogram_drawing(MeticsValues1,  outputPathname, outputBasicMetricschartPath);
		Generatehtml(outputBasicMetricschart);
		cout << endl << "Details of the assessment results are saved in:" << outputPathname + htmlFilename<<endl;
	}
}

//Compare multiple data sets
void multipledataset(vector< vector<float> > MeticsValues, vector<string> &sv_files1, vector<string> &tool_names, string &outputBasicMetricschart, string &newInfo){
    // Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
	vector<string> scenarios;
	if(tool_names.size()>1) scenarios = tool_names;
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			replaceUnderscoreWithDot(str);
			scenarios.push_back(str);
		}
	}
    // Create a data file to hold the data
	string filename = "evaluation_result";
	string filenamePath = outputBasicMetricschart + '/'+ filename;
    ofstream dataFile(filenamePath);
    if (!dataFile.is_open()) {
        cerr << "Unable to create data file: " << filenamePath << endl;
        return ;
    }

    // Write data to a file
    for (size_t i = 0; i < scenarios.size(); i++) {
        dataFile << scenarios[i];
        for (size_t j = 0; j < 3; j++) {
            dataFile << " " << MeticsValues[i][j];
        }
        dataFile << endl;
    }
    dataFile.close();

    // Set the file name to save
    string outputFileName = "evaluation_result.png";
    string outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;

    // Create a pipe connection to GNUplot
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        cerr << "Unable to open the GNUplot pipeline" << endl;
        return ;
    }

    // Use the GNUplot command to plot the bar plot and save the output as a PNG file
    fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
    fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800\n");
    fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
    fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
    fprintf(gnuplotPipe, "set title 'Performance comparison between different tools' font 'Times-Roman,14'\n");
    fprintf(gnuplotPipe, "set xlabel 'Tool name' font 'Times-Roman,12'\n");
    fprintf(gnuplotPipe, "set ylabel 'Percentage' font 'Times-Roman,12'\n");
    fprintf(gnuplotPipe, "set style fill solid\n");
    fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
    fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,10'\n");
    fprintf(gnuplotPipe, "set yrange [0:1]\n");
    fprintf(gnuplotPipe, "set ytics 0.1\n");
    fprintf(gnuplotPipe, "set style data histograms\n");

    // Draw a bar chart, specifying the number of columns for each data column purple  olive magenta
    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb '#568AC6' title 'Recall', '' using 3 lc rgb '#F08700' title 'Precision', '' using 4 lc rgb '#F06767' title 'F1'\n", filenamePath.c_str());

    // Close the GNUplot pipeline
    pclose(gnuplotPipe);

    AddfileInformation(filenamePath, newInfo);
    cout << "Performance metrics results between different data sets are saved in: " << outputFileNamePath << endl;

}

//stat 5
void multipledataset(vector< vector<int> > MeticsValues, vector<string> &sv_files1, vector<string> &tool_names, string &outputBasicMetricschart, string &newInfo){
    // Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
	vector<string> scenarios;
	if(tool_names.size()>1) scenarios = tool_names;
	//Note to be modified to remove '/'
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			replaceUnderscoreWithDot(str);
			scenarios.push_back(str);
		}
	}
    // Create a data file to hold the data
	string filename = "result_classification";
	string filenamePath = outputBasicMetricschart + '/'+ filename;
    ofstream dataFile(filenamePath);
    if (!dataFile.is_open()) {
        cerr << "Unable to create data file: " << filenamePath << endl;
        return ;
    }
    // Write data to a file
    for (size_t i = 0; i < scenarios.size(); i++) {
        dataFile << scenarios[i];
        for (size_t j = 0; j < 4; j++) {
            dataFile << " " << MeticsValues[i][j];
        }
        dataFile << endl;
    }
    dataFile.close();

    // Set the file name to save
    string outputFileName = "result_classification.png";
    string outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;

    // Create a pipe connection to GNUplot
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        cerr << "Unable to open the GNUplot pipeline" << endl;
        return ;
    }

    // Use the GNUplot command to plot the bar plot and save the output as a PNG file
    fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
    fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800\n");
    fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
    fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
    fprintf(gnuplotPipe, "set title 'Benchmark results between different tools' font 'Times-Roman,14'\n");
    fprintf(gnuplotPipe, "set xlabel 'Tool name' font 'Times-Roman,12'\n");
    fprintf(gnuplotPipe, "set ylabel 'Count' font 'Times-Roman,12'\n");
    fprintf(gnuplotPipe, "set style fill solid\n");
    fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
    fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,10'\n");
    fprintf(gnuplotPipe, "set yrange [0:*]\n");
    fprintf(gnuplotPipe, "set style data histograms\n");

    // Draw a bar chart, specifying the number of columns for each data column purple blue green red
    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb '#6FBE47' title 'TP.user', '' using 3 lc rgb '#FF6847' title 'TP.bench', '' using 4 lc rgb '#FDD700' title 'FP', '' using 5 lc rgb '#AB6520' title 'FN'\n", filenamePath.c_str());
//    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'blue' title 'Recall', '' using 3 lc rgb 'green' title 'Precision', '' using 4 lc rgb 'red' title 'F1'\n", outputFileNamePath.c_str());

    // Close the GNUplot pipeline
    pclose(gnuplotPipe);
    AddfileInformation(filenamePath, newInfo);
    cout << endl << "The classification of results between different data sets is saved in: " << outputFileNamePath << endl;
}

//Modifying File Information
void AddfileInformation(string &FileNamePath, string &Info){
	ifstream inputFile(FileNamePath);
	if (!inputFile.is_open()) {
		cerr << "Unable to open file" << endl;
		return ;
	}
	// Read the original contents of the file
	vector<string> fileContents;
	string line;
	while (getline(inputFile, line)) {
		fileContents.push_back(line);
	}
	inputFile.close();
	string tempFilename = FileNamePath + ".temp";
	ofstream tempFile(tempFilename);
	if (!tempFile.is_open()) {
		std::cerr << "Unable to create temporary file" << std::endl;
		return ;
	}
	// Writes a new header line
	tempFile << Info << endl;
	// Write the contents of the original file back to the temporary file, keeping the original contents
	for (const std::string& content : fileContents) {
		tempFile << content << std::endl;
	}
	tempFile.close();
	// Delete original file
	if (remove(FileNamePath.c_str()) != 0) {
		std::cerr << "Unable to delete the original file" << std::endl;
		return ;
	}
	// Rename the temporary file to the original file
	if (rename(tempFilename.c_str(), FileNamePath.c_str()) != 0) {
		std::cerr << "Cannot rename temporary files to original files" << std::endl;
		return ;
	}
}

//stat 2
void Histogram_drawing(vector< vector<float> > MeticsValues, string &outputPathname, string &outputBasicMetricschart){
	vector<pair<string, float>> dataV;
	dataV.push_back(make_pair("Recall", MeticsValues[0][0]));
	dataV.push_back(make_pair("F1-score", MeticsValues[0][1]));
	dataV.push_back(make_pair("Precision", MeticsValues[0][2]));

	string outputFileName = "evaluation_result.png";
	string outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;

	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}
	fprintf(gnuplotPipe, "set term pngcairo size 800,600\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title 'Performance Metric Results Statistics' font 'Times-Roman,16'\n");
	fprintf(gnuplotPipe, "set xlabel 'performance metric' font 'Times-Roman,14'\n");
	fprintf(gnuplotPipe, "set ylabel 'Percentage'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.2\n");
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");
	fprintf(gnuplotPipe, "plot '-' using 2:xtic(1) title 'Percentage'\n");

	for (const auto& pair : dataV) {
		fprintf(gnuplotPipe, "%s %lf\n", pair.first.c_str(), pair.second);
	}
	fprintf(gnuplotPipe, "e\n");
	pclose(gnuplotPipe);
	cout << "Performance metrics results between different data sets are saved in: " << outputFileNamePath << endl;
}

// stat 2
void Histogram_drawing(vector< vector<int> > MeticsValues1, string &outputPathname, string &outputBasicMetricschart){
		string outputFileName, outputFileNamePath;
		// Prepare the bar chart data
	    vector<pair<string, float>> dataV;
	    dataV.push_back(make_pair("TP.user", MeticsValues1[0][0]));
	    dataV.push_back(make_pair("TP.bench", MeticsValues1[0][1]));
	    dataV.push_back(make_pair("FP", MeticsValues1[0][2]));
	    dataV.push_back(make_pair("FN", MeticsValues1[0][3]));

	    outputFileName = "result_classification.png";
	    outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;
	    // Create a pipe connection to GNUplot
	    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	    if (!gnuplotPipe) {
	        cerr << "Unable to open the GNUplot pipeline" << endl;
	        return ;
	    }
	    // Use the GNUplot command to plot the bar plot and save the output as a PNG file
	    fprintf(gnuplotPipe, "set term pngcairo size 800,600\n");
	    fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str());
	    fprintf(gnuplotPipe, "set title 'Statistics of Benchmark Results' font 'Times-Roman,16'\n");
	    fprintf(gnuplotPipe, "set xlabel 'benchmark evaluation metric' font 'Times-Roman,14'\n");
	    fprintf(gnuplotPipe, "set ylabel 'Count'\n");
	    fprintf(gnuplotPipe, "set style fill solid\n");
	    fprintf(gnuplotPipe, "set boxwidth 1.5\n");
	    fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
	    fprintf(gnuplotPipe, "set yrange [0:*]\n");
	    fprintf(gnuplotPipe, "set style data histograms\n");
	    fprintf(gnuplotPipe, "plot '-' using 2:xtic(1) title 'count'\n");
	    // Pass the data to GNUplot
	    for (const auto& pair : dataV) {
	        fprintf(gnuplotPipe, "%s %lf\n", pair.first.c_str(), pair.second);
	    }
	    fprintf(gnuplotPipe, "e\n");
	    pclose(gnuplotPipe);
	    cout << endl << "The classification of results between different data sets is saved in " << outputFileNamePath << endl;
}

// stat 1
void SvNumberDistributionGraph(int max_valid_reg_thres, string &refRegSizeFinename, string &refRegSizeFinename_tmp){
    string fileName, fileNamePath, filePathstr_tmp;
    if(refRegSizeFinename.compare("ref_reg_size_benchmark")==0){
	    fileName = refRegSizeFinename;
	    fileNamePath = refRegSizeFinename_tmp;
	    filePathstr_tmp = getContentAfterSlash(fileNamePath) + ".png";
	    if(folderPng1.size()==0) folderPng1.push_back(filePathstr_tmp);
    }else if(refRegSizeFinename.compare("ref_reg_size_user")== 0 and max_valid_reg_thres == 0){
	    fileName = refRegSizeFinename;
	    fileNamePath = refRegSizeFinename_tmp;
	    filePathstr_tmp = getContentAfterSlash(fileNamePath) + ".png";
	    folderPng1.push_back(filePathstr_tmp);
    }else if(refRegSizeFinename.compare("ref_reg_size_user")== 0 and max_valid_reg_thres > 0){
	    fileName = refRegSizeFinename + "_long_filtered";
	    fileNamePath = refRegSizeFinename_tmp;
	    filePathstr_tmp = getContentAfterSlash(fileNamePath) + ".png";
	    folderPng1.push_back(filePathstr_tmp);
    }

	ifstream inputFile(fileNamePath.c_str(), ios::in);
	if (!inputFile) {
		cerr << "Unable to open data file: " << fileNamePath << endl;
		return ;
	}
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");

	if (!gnuplotPipe) {
		cerr << "Unable to start GNUplot." << endl;
		return ;
	}
	// Set the GNUplot option
	fprintf(gnuplotPipe, "set title 'Quantitative distribution'\n");
	fprintf(gnuplotPipe, "set xlabel 'SV size'\n");
	fprintf(gnuplotPipe, "set ylabel 'Count'\n");
	fprintf(gnuplotPipe, "set xrange [-100:100]\n");
	fprintf(gnuplotPipe, "set terminal png  size 800,600\n");
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set output '%s.png'\n", fileNamePath.c_str());

	// Draw the graph with the filename variable as the data filename
	fprintf(gnuplotPipe, "plot '%s' every ::2 using 1:2 with linespoints pt 2 linewidth 2 title 'count'\n", fileNamePath.c_str());
	// Close GNUplot
	fflush(gnuplotPipe);
	pclose(gnuplotPipe);

	inputFile.close();
}

//stat 3
void CenterdistanceAndAreasizeratio(string &sizeDifStatDirname){
	string filename1 = "sv_size_ratio_stat_long_filtered", filename2 = "sv_size_dif_stat_long_filtered", filename_pngPath_tmp;
	string title1 = "SV size difference ratio", title2 = "Center distance difference";
	string fileNamePath1 = sizeDifStatDirname + filename1, fileNamePath2 = sizeDifStatDirname + filename2;
	string filename_png = "result_details.png";
	string filename_pngPath = sizeDifStatDirname + filename_png;
	SVSizeRatioFile.push_back(fileNamePath1);
	filename_pngPath_tmp = getContentAfterSlash(filename_pngPath);
	folderPng3.push_back(filename_pngPath_tmp);
	// Create a pipe connection to the GNUplot process
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline." << endl;
		return ;
	}
	// Read the first data file
	ifstream dataFile1(fileNamePath1);

	if (!dataFile1.is_open()) {
		cerr << "Unable to open the first data file: " << fileNamePath1 << endl;
		return ;
	}

	// Prepare a container for storing data
	vector<string> xLabels1;
	vector<double> yValues1;
	// Read data
	string line1;
	while (getline(dataFile1, line1)) {
		istringstream iss(line1);
		string xLabel;
		double yValue;

		if (iss >> xLabel >> yValue) {
			xLabels1.push_back(xLabel);
			yValues1.push_back(yValue);
		}
	}
	dataFile1.close();

	ifstream dataFile2(fileNamePath2);
	if (!dataFile2.is_open()) {
		cerr << "Unable to open second data file: " << fileNamePath2 << endl;
		return ;
	}

	vector<string> xLabels2;
	vector<double> yValues2;

	string line2;
	while (getline(dataFile2, line2)) {
		istringstream iss(line2);
		string xLabel;
		double yValue;

		if (iss >> xLabel >> yValue) {
			xLabels2.push_back(xLabel);
			yValues2.push_back(yValue);
		}
	}
	dataFile2.close();
	stringstream gnuplotCommands;
	gnuplotCommands << "set terminal pngcairo size 1500,1200\n";
	gnuplotCommands << "set title '" << title1.c_str() << "' font 'Times-Roman,20'\n";
	gnuplotCommands << "set output '" << filename_pngPath.c_str()<< "'\n";
	gnuplotCommands << "set style data histograms\n";
	gnuplotCommands << "set style fill solid border -1\n";
	gnuplotCommands << "set border linecolor rgb 'black' linewidth 2\n";
	gnuplotCommands << "set boxwidth 0.5 relative\n";
	gnuplotCommands << "set xtics nomirror\n";
	gnuplotCommands << "set yrange [0:*]\n";
	gnuplotCommands << "set multiplot layout 2, 1\n";
	gnuplotCommands << "plot '-' using 2:xtic(1) with boxes lc rgb '#2278B4' title 'count'\n";

	fprintf(gnuplotPipe, "%s", gnuplotCommands.str().c_str());
	fflush(gnuplotPipe);
	// Send the data point of the first data file to GNUplot
	for (size_t i = 0; i < xLabels1.size(); ++i) {
		fprintf(gnuplotPipe, "'%s' %f\n", xLabels1[i].c_str(), yValues1[i]);
	}
	fprintf(gnuplotPipe, "e\n"); // End the first data point input

	// Plot the second histogram in multiplot
	gnuplotCommands.str(""); // Clear previous commands
	gnuplotCommands << "set title '" << title2.c_str() << "' font 'Times-Roman,20'\n";
	gnuplotCommands << "set xtics nomirror rotate by -45\n";
	gnuplotCommands << "set boxwidth 1.0 relative\n"; // Set the column width for the second diagram
	gnuplotCommands << "set border linecolor rgb 'black' linewidth 2\n";
//    gnuplotCommands << "set label 'Center distance difference' at graph 0.3, 0.9 left font 'Time New Roman,18\n"; // Add diagram notes to the second diagram
	gnuplotCommands << "plot '-' using 2:xtic(1) with boxes lc rgb '#E20613'  title 'count'\n";
	fprintf(gnuplotPipe, "%s", gnuplotCommands.str().c_str());
	fflush(gnuplotPipe);

	// Send the data points of the second data file to GNUplot
	for (size_t i = 0; i < xLabels2.size(); ++i) {
		fprintf(gnuplotPipe, "'%s' %f\n", xLabels2[i].c_str(), yValues2[i]);
	}
	fprintf(gnuplotPipe, "e\n"); // End the second data point input

	pclose(gnuplotPipe);

	cout << "\n" <<"Statistical results of center distance and area size ratio are saved in:" << filename_pngPath << endl;
}

//stat 4
void SVsizeAndNumstatistics(string &sizeNumStatDirname, vector< vector<float> > MeticsValues_4){
	string filename, filenamePath, outputFileName, outputFileNamePath, outputFileNamePath_tmp;
	// Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
	vector<string> scenarios ={"1-100bp","101-250bp","251-500bp","501-1000bp","1001-5000bp","5001-10000bp",">10000bp"};

	// Create a data file to hold the data
	filename = "calculation_result";
	filenamePath = sizeNumStatDirname + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file: " << filenamePath << endl;
		return ;
	}

	// Write data to a file
	for (size_t i = 0; i < scenarios.size(); i++) {
		dataFile << scenarios[i];
		for (size_t j = 0; j < 3; j++) {
			dataFile << " " << MeticsValues_4[i][j];
		}
		dataFile << endl;
	}
	dataFile.close();

	// Set the file name to save
	outputFileName = "calculation_result.png";
	outputFileNamePath = sizeNumStatDirname + '/' + outputFileName;
	outputFileNamePath_tmp = getContentAfterSlash(outputFileNamePath);
	folderPng4.push_back(outputFileNamePath_tmp);
	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}

	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title 'Number statistics of classified metrics in different SV ranges' font 'Times-Roman,14'\n");
	fprintf(gnuplotPipe, "set xlabel 'SV ranges'\n");
	fprintf(gnuplotPipe, "set ylabel 'Percentage'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,10'\n");
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");

	// Draw a bar chart, specifying the number of columns for each data column
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb '#3E4998' title 'Recall', '' using 3 lc rgb '#D95F02' title 'Precision', '' using 4 lc rgb '#E7298A' title 'F1'\n", filenamePath.c_str());
//    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'blue' title 'Recall', '' using 3 lc rgb 'green' title 'Precision', '' using 4 lc rgb 'red' title 'F1'\n", outputFileNamePath.c_str());
	// Close the GNUplot pipeline
	pclose(gnuplotPipe);
	string newInfo = "Range Recall Precision F1-score";
	AddfileInformation(filenamePath, newInfo);
	cout << endl << "Statistical results of indicators with different SV sizes are saved in: " << outputFileNamePath << endl;
}

//stat 4
void SVsizeAndNumstatistics(string &sizeNumStatDirname, vector< vector<int> > MeticsValues1_4){
	string filename, filenamePath, outputFileName, outputFileNamePath, outputFileNamePath_tmp;
// Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
	vector<string> scenarios ={"1-100bp","101-250bp","251-500bp","501-1000bp","1001-5000bp","5001-10000bp",">10000bp"};

	// Create a data file to hold the data
	filename = "quantitative_statistics";
	filenamePath = sizeNumStatDirname + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file: " << filenamePath << endl;
		return ;
	}

	// Write data to a file
	for (size_t i = 0; i < scenarios.size(); i++) {
		dataFile << scenarios[i];
		for (size_t j = 0; j < 4; j++) {
			dataFile << " " << MeticsValues1_4[i][j];
		}
		dataFile << endl;
	}
	dataFile.close();

	// Set the file name to save
	outputFileName = "quantitative_Statistics.png";
	outputFileNamePath = sizeNumStatDirname + '/' + outputFileName;
	outputFileNamePath_tmp=getContentAfterSlash(outputFileNamePath);
	folderPng4.push_back(outputFileNamePath_tmp);
	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}

	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title 'Quantitative Statistics'\n");
	fprintf(gnuplotPipe, "set xlabel 'SV Ranges'\n");
	fprintf(gnuplotPipe, "set ylabel 'Count'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,10'\n");
	fprintf(gnuplotPipe, "set yrange [0:*]\n");
	fprintf(gnuplotPipe, "set style data histograms\n");

	// Draw a bar chart, specifying the number of columns for each data column
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb '#A09DCB' title 'TP.user', '' using 3 lc rgb '#73B8FA' title 'TP.bench', '' using 4 lc rgb '#69BB1E' title 'FP', '' using 5 lc rgb '#EB6365' title 'FN'\n", filenamePath.c_str());
//	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'olive' title 'Precision', '' using 4 lc rgb 'magenta' title 'F1'\n", filenamePath.c_str());

	// Close the GNUplot pipeline
	pclose(gnuplotPipe);
	string newInfo = "Range TP_user TP_bench FP FN";
	AddfileInformation(filenamePath, newInfo);
	cout << "Quantitative statistics of indicators with different SV sizes are saved in: " << outputFileNamePath << endl << endl;
}

void ComparisonofMetricsindifferentranges(string &DiffRangeStatDirname, vector< vector<float> > meticsvalues, vector<string> &tool_names, string &FileName, string &outputfilename, vector<string> &sv_files1){
// Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
	vector<string> scenarios;
	if(tool_names.size()>1) scenarios = tool_names;
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			scenarios.push_back(str);
		}
	}
	string title = "Range:" + FileName;

	// Create a data file to hold the data
	string filename = FileName;
	string filenamePath = DiffRangeStatDirname + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file" << endl;
		return ;
	}

	// Write data to a file
	for (size_t i = 0; i < scenarios.size(); i++) {
		dataFile << scenarios[i];
		for (size_t j = 0; j < 3; j++) {
			dataFile << " " << meticsvalues[i][j];
		}
		dataFile << endl;
	}
	dataFile.close();

	// Set the file name to save
	string outputFileName = outputfilename;
	string outputFileNamePath = DiffRangeStatDirname + '/' + outputFileName;

	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}

	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 800,600\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", title.c_str());
	fprintf(gnuplotPipe, "set xlabel 'Tool Name'\n");
	fprintf(gnuplotPipe, "set ylabel 'Value'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
//    fprintf(gnuplotPipe, "set xrange [-0.5:2.5]\n"); // Change the X-axis range
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");

	// Draw a bar chart, specifying the number of columns for each data column
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", filenamePath.c_str());
//    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'blue' title 'Recall', '' using 3 lc rgb 'green' title 'Precision', '' using 4 lc rgb 'red' title 'F1'\n", outputFileNamePath.c_str());
	// Close the GNUplot pipeline
	pclose(gnuplotPipe);
	string newInfo = "tool Recall Precision F1-score";
	AddfileInformation(filenamePath, newInfo);
	cout << endl << "The comparison results of different SV range indicators are saved in: " << outputFileNamePath << endl;
}

void ComparisonofMetricsindifferentrangesOp(string &outputDiffRangeBasicMetricschart, vector<string> &tool_names, vector<string> &sv_files1){
	string range = "1-100",range_0 = "Range1-100.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_0, tool_names, range, range_0, sv_files1);
	string range1 = "101-250",range_1 = "Range101-250.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_1, tool_names, range1, range_1, sv_files1);
	string range2 = "251-500",range_2 = "Range251-500.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_2, tool_names, range2, range_2, sv_files1);
	string range3 = "501-1000",range_3 = "Range501-1000.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_3, tool_names, range3, range_3, sv_files1);
	string range4 = "1001-5000",range_4 = "Range1001-5000.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_4, tool_names, range4, range_4, sv_files1);
	string range5 = "5001-10000",range_5 = "Range5001-10000.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_5, tool_names, range5, range_5, sv_files1);
	string range6 = "GreaterThan10001",range_6 = "GreaterThan10001.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_6, tool_names, range5, range_6, sv_files1);
}

//stat 5
void GenerateFileOp(string &DiffRangeStatDirname, vector<string> &tool_names, vector<string> &filenames, vector<string> &sv_files1){
	string FileNameStr0 = "1-100";
	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr0, MeticsValues4_0, filenames, sv_files1);
	string FileNameStr1 = "101-250";
	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr1, MeticsValues4_1, filenames, sv_files1);
	string FileNameStr2 = "251-500";
	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr2, MeticsValues4_2, filenames, sv_files1);
	string FileNameStr3 = "501-1000";
	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr3, MeticsValues4_3, filenames, sv_files1);
	string FileNameStr4 = "1001-5000";
	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr4, MeticsValues4_4, filenames, sv_files1);
	string FileNameStr5 = "5001-10000";
	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr5, MeticsValues4_5, filenames, sv_files1);
	string FileNameStr6 = "greaterthan10001";
	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr6, MeticsValues4_6, filenames, sv_files1);
}
//stat 5
void GenerateFile(string &DiffRangeStatDirname, vector<string> &tool_names, string &FileName, vector< vector<float> > meticsvalues, vector<string> &filenames, vector<string> &sv_files1){
	vector<string> scenarios;
	if(tool_names.size()>1) scenarios = tool_names;
	//Note to be modified to remove '/'’
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			replaceUnderscoreWithDot(str);
			scenarios.push_back(str);
		}
	}
	string title = "Range:" + FileName;
	// Create a data file to hold the data
	string filename = FileName;
	string filenamePath = DiffRangeStatDirname + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file: " << filenamePath << endl;
		return ;
	}
	// Write data to a file
	for (size_t i = 0; i < scenarios.size(); i++) {
		dataFile << scenarios[i];
		for (size_t j = 0; j < 3; j++) {
			dataFile << " " << meticsvalues[i][j];
		}
		dataFile << endl;
	}
	dataFile.close();
	filenames.push_back(filenamePath);
}

//NEW stat 5
void GenerateMultiBarCharts(string &outputBasicMetricschart, vector<string>& fileNames, string Header_line) {
	string outputFileName, outputFileNamePath, outputFileNamePath_tmp;
	outputFileName = "different_range.png";
	outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;
	outputFileNamePath_tmp=getContentAfterSlash(outputFileNamePath);
	folderPng5.push_back(outputFileNamePath_tmp);
	vector<string> titles =  {
			        "1-100 bp", "101-250 bp", "251-500 bp",
			        "501-1000 bp", "1001-5000 bp", "5001-10000 bp", ">=10001 bp"
			    };
	//// Open the Gnuplot pipeline
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        cerr << "Unable to start GNUplot." << endl;
        return;
    }
	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
    //Fig1
	fprintf(gnuplotPipe, "set term pngcairo size 1500,1800\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set multiplot layout 4, 2\n");

	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[0].c_str());
	fprintf(gnuplotPipe, "set xlabel 'Tool Name'\n");
	fprintf(gnuplotPipe, "set ylabel 'Percentage'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");

	// Draw a bar chart, specifying the number of columns for each data column
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", fileNames[0].c_str());
//    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'blue' title 'Recall', '' using 3 lc rgb 'green' title 'Precision', '' using 4 lc rgb 'red' title 'F1'\n", outputFileNamePath.c_str());

	//Fig2
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[1].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", fileNames[1].c_str());

	//Fig3
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[2].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", fileNames[2].c_str());

	//Fig4
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[3].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", fileNames[3].c_str());

	//Fig5
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[4].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", fileNames[4].c_str());

	//Fig6
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[5].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", fileNames[5].c_str());

	//Fig7
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[6].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", fileNames[6].c_str());

	// Disable the multi-image mode
    fprintf(gnuplotPipe, "unset multiplot\n");

    // Close Gnuplot
    fflush(gnuplotPipe);
    pclose(gnuplotPipe);
    cout << endl << "The comparison of metrics in different ranges is saved in:" <<outputFileNamePath<< endl;
    for(size_t i=0; i<fileNames.size(); i++){
    	AddfileInformation(fileNames[i], Header_line);
    }
}

//stat 5
//void GenerateSVsizeRatioFileGraph(vector<string> svfilesV, vector<string> toolnameV, string &FileSavingPath){
//	 vector<string> titles;
//	 string filename = "sv_size_ratio.png", filenamePath;
//	 filenamePath = FileSavingPath + '/' + filename;
//	 if(toolnameV.size()>0)	titles = toolnameV;
//	 else {
//	 		for(size_t i=0; i<svfilesV.size();i++){
//	 			size_t lastSlashPos = svfilesV.at(i).find_last_of('/');
//	 			string str;
//	 			if (lastSlashPos == string::npos) {
//	 				// If '/' is not found, the entire string is printed
//	 				 str = svfilesV.at(i-1);
//	 			} else {
//	 				// When '/' is reached, the part after '/' is printed
//	 				str = svfilesV.at(i).substr(lastSlashPos + 1);
//	 			}
//	 			titles.push_back(str);
//	 		}
//	 	}
////	 else titles = svfilesV;
//	 // Create a GNU Plot pipeline
//	 FILE* gnuplotPipe = popen("gnuplot -persist", "w");
//	 if (!gnuplotPipe) {
//		 std::cerr << "Unable to start GNU Plot." << std::endl;
//		 return ;
//	 }
//
//	 fprintf(gnuplotPipe, "set terminal pngcairo size 1200,1000\n");
//	 fprintf(gnuplotPipe, "set output '%s'\n", filenamePath.c_str());
//
//	 // Enter multiplot mode and set the layout
//	 fprintf(gnuplotPipe, "set multiplot layout %d,1\n", (int)SVSizeRatioFile.size());
//
//	 // Loop the subgraph
//	 for (size_t i = 0; i < SVSizeRatioFile.size(); ++i) {
//		 fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[i].c_str());
//		 fprintf(gnuplotPipe, "set style data histograms\n");
//		 fprintf(gnuplotPipe, "set style fill solid border -1\n");
//		 fprintf(gnuplotPipe, "set boxwidth 0.5 relative\n");
//		 fprintf(gnuplotPipe, "set xtics nomirror\n");
//		 fprintf(gnuplotPipe, "set yrange [0:*]\n");
//		 fprintf(gnuplotPipe, "set xlabel 'SV size ratio range'\n");
//		 fprintf(gnuplotPipe, "set ylabel 'Count'\n");
//
//		 // subgraph
//		 fprintf(gnuplotPipe, "plot '-' using 2:xtic(1) with boxes title 'count'\n");
//
//		 ifstream dataFile(SVSizeRatioFile[i]);
//		 if (dataFile.is_open()) {
//			 string line;
//			 while (getline(dataFile, line)) {
//				 fprintf(gnuplotPipe, "%s\n", line.c_str());
//			 }
//			 dataFile.close();
//		 }
//		 fprintf(gnuplotPipe, "e\n"); // The subgraph is drawn
//	 }
//
//	 //Exit the multiplot mode
//	 fprintf(gnuplotPipe, "unset multiplot\n");
//	 // Close the GNU Plot pipeline
//	 pclose(gnuplotPipe);
//
//	 //Print file path
//	 cout << endl << "The results of the sv size ratio identified by different tools are saved in: " << filenamePath << endl;
//}

//Replace above function: stat 5
void GenerateSVsizeRatioFileGraph(vector<string> svfilesV, vector<string> toolnameV, string &FileSavingPath) {
    vector<string> annotations;
    string filename = "sv_size_ratio", filename1 = "sv_quantity_ratio", filenamePath, filename1Path, outputFileNamePath, outputFileName = "sv_size_ratio.png";
    filenamePath = FileSavingPath + '/' + filename;
    filename1Path =  FileSavingPath + '/' + filename1;
    outputFileNamePath = FileSavingPath + '/' +outputFileName;

    if (toolnameV.size() > 0) {
    	annotations = toolnameV;
    } else {
        for (size_t i = 0; i < svfilesV.size(); i++) {
            size_t lastSlashPos = svfilesV.at(i).find_last_of('/');
            string str;
            if (lastSlashPos == string::npos) {
                str = svfilesV.at(i - 1);
            } else {
                str = svfilesV.at(i).substr(lastSlashPos + 1);
            }
            replaceUnderscoreWithDot(str);
            annotations.push_back(str);
        }
    }
   vector<string> header_line ={"0.0-0.5","0.5-0.7","0.7-1.2","1.2-2.0","2.0-5.0","5.0-10.0","10.0-50.0","50.0-100.0",">100.0"};
   ofstream outputFile(filenamePath);

   if (!outputFile.is_open()) {
	   cerr << "Unable to open file: " << filenamePath << endl;
	   return ;
   }
   // Iterate over containers and write their values and header information to a file
   for (size_t i = 0; i < SizeRatioV.size(); ++i) {
	   outputFile << header_line[i] << " "; // Add title information

	   for (int value : SizeRatioV[i]) {
		   outputFile << value << " ";
	   }
	   outputFile << endl; // Insert newlines between containers
   }
   outputFile.close();
   cout << endl << "The sv size ratio of the text information of different tools is saved to: " << filenamePath << endl;

   ofstream outputFile1(filename1Path);

  if (!outputFile1.is_open()) {
   cerr << "Unable to open file: " << filename1Path << endl;
   return ;
  }
  // Iterate over containers and write their values and header information to a file
  for (size_t i = 0; i < SizeRatio_V.size(); ++i) {
   outputFile1 << header_line[i] << " "; // Add title information

   for (double value : SizeRatio_V[i]) {
	   outputFile1 << value << " ";
   }
   outputFile1 << endl; // Insert newlines between containers
  }
  outputFile1.close();
  cout << endl << "The ratio of overlapping SVs of different sizes identified by different tools is saved to: " << filename1Path << endl;

   // Create a pipe connection to GNUplot
   FILE* gnuplotPipe = popen("gnuplot -persist", "w");
   if (!gnuplotPipe) {
	   cerr << "Unable to open the GNUplot pipeline" << endl;
	   return ;
   }
   fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
   fprintf(gnuplotPipe, "set terminal pngcairo size 1200,800\n");
   fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
   fprintf(gnuplotPipe, "set title 'Quantity Statistics of Overlapping SV Size Ratios from Different Tools' font 'Times-Roman,14'\n");
   fprintf(gnuplotPipe, "set xlabel 'Range of size ratio'\n");
   fprintf(gnuplotPipe, "set ylabel 'Count (histogram)'\n");
   fprintf(gnuplotPipe, "set style fill solid\n");
   fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
   fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
   fprintf(gnuplotPipe, "set ytics nomirror\n");
   fprintf(gnuplotPipe, "set yrange [0:*]\n");

   fprintf(gnuplotPipe, "set y2tics\n");
   fprintf(gnuplotPipe, "set y2label 'Percentage (line)'\n");
   fprintf(gnuplotPipe, "set y2range [0:1]\n");
   fprintf(gnuplotPipe, "set y2tics 0.1\n");

   // Style the bar chart based on the number of input data sets
   if(svfilesV.size()==2){
//	   fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title '%s', '' using 3 lc rgb 'red' title '%s' \n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str());
	   fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str());
	   fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 \n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str());
	   cout << endl << "The identification results of the SVs region size ratio of different tools are saved in" << outputFileNamePath << endl;
   }else if(svfilesV.size()==3){
	   fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str());
	   fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str());
	   cout << endl << "The identification results of the SVs region size ratio of different tools are saved in" << outputFileNamePath << endl;
   }else if(svfilesV.size()==4){
	   fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str());
	   fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str());
	   cout << endl << "The identification results of the SVs region size ratio of different tools are saved in" << outputFileNamePath << endl;
   }else if(svfilesV.size()==5){
	   fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram, '' using ($6+0.2) lc rgb '#FAB833' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str());
	   fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2, '' using 6 lc rgb '#FAB833' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str());
	   cout << endl << "The identification results of the SVs region size ratio of different tools are saved in" << outputFileNamePath << endl;
   }else if(svfilesV.size()==6){
	   fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram, '' using ($6+0.2) lc rgb '#FAB833' title '%s' with histogram, '' using ($7+0.2) lc rgb '#227CDE' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str());
	   fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2, '' using 6 lc rgb '#FAB833' lw 2 title '%s' with linespoints axes x1y2, '' using 7 lc rgb '#227CDE' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str());
	   cout << endl << "The identification results of the SVs region size ratio of different tools are saved in" << outputFileNamePath << endl;
   }else if(svfilesV.size()==7){
	   fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram, '' using ($6+0.2) lc rgb '#FAB833' title '%s' with histogram, '' using ($7+0.2) lc rgb '#227CDE' title '%s' with histogram, '' using ($8+0.2) lc rgb '#55436E' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str(), annotations[6].c_str());
	   fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2, '' using 6 lc rgb '#FAB833' lw 2 title '%s' with linespoints axes x1y2, '' using 7 lc rgb '#227CDE' lw 2 title '%s' with linespoints axes x1y2, '' using 8 lc rgb '#55436E' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str(), annotations[6].c_str());
	   cout << endl << "The identification results of the SVs region size ratio of different tools are saved in" << outputFileNamePath << endl;
   }else if(svfilesV.size()==8){
	   fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram, '' using ($6+0.2) lc rgb '#FAB833' title '%s' with histogram, '' using ($7+0.2) lc rgb '#227CDE' title '%s' with histogram, '' using ($8+0.2) lc rgb '#55436E' title '%s' with histogram, '' using ($9+0.2) lc rgb '#A7C4E1' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str(), annotations[6].c_str(), annotations[7].c_str());
	   fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2, '' using 6 lc rgb '#FAB833' lw 2 title '%s' with linespoints axes x1y2, '' using 7 lc rgb '#227CDE' lw 2 title '%s' with linespoints axes x1y2, '' using 8 lc rgb '#55436E' lw 2 title '%s' with linespoints axes x1y2, '' using 9 lc rgb '#A7C4E1' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str(), annotations[6].c_str(), annotations[7].c_str());
	   cout << endl << "The identification results of the SVs region size ratio of different tools are saved in" << outputFileNamePath << endl;
   }else
	   cout << endl << "Up to eight evaluation tools for region size ratio images, exceeding six won't create images." << endl;

   string Add_Info = "#sizeRatio";
   for(size_t i=0; i<annotations.size(); i++){
	   Add_Info = Add_Info + " " + annotations[i];
   }
   AddfileInformation(filenamePath, Add_Info);
   AddfileInformation(filename1Path, Add_Info);
}

