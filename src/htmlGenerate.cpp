#include "htmlGenerate.h"

void Generatehtml(string figuresFilePath){
	int32_t table_num = 1, figure_num = 1;
	string HtmlFilePath = outputPathname + '/' + htmlFilename;
	ofstream htmlFile(HtmlFilePath);

	if (htmlFile.is_open()) {
		htmlFile << "<!DOCTYPE html>\n";
		htmlFile << "<html>\n";
		htmlFile << "<head>\n";
		htmlFile << "  <title>Evaluation Reports</title>\n";
		htmlFile << "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
		htmlFile << "  <title>Table Example</title>\n";
		htmlFile << "  <style>\n";
		// Add CSS styles for titles
		htmlFile << "    body {\n";
		htmlFile << "      font-family: 'Arial', sans-serif; /* Set body font family */\n";
		htmlFile << "      line-height: 1.6; /* Set line height for better readability */\n";
		htmlFile << "      margin: 0; /* Remove default margin */\n";
		htmlFile << "      padding: 0; /* Remove default padding */\n";
		htmlFile << "      background-color: #f9f9f9; /* Set background color */\n";
		htmlFile << "      color: #333; /* Set text color */\n";
		htmlFile << "    }\n";
		htmlFile << "    pre{\n";
		htmlFile << "      background - color: #f4f4f4;\n";
		htmlFile << "      padding: 10px;\n";
		htmlFile << "      border: 1px solid #ccc;\n";
		htmlFile << "      border: 1px solid #ccc;\n";
		htmlFile << "      border - radius: 5px;\n";
		htmlFile << "      white - space: pre - wrap;\n";
		htmlFile << "      font - family: 'Courier New', monospace;\n";
		htmlFile << "      font-size: 16px;\n";
		htmlFile << "    }\n";
		htmlFile << "   .command-container {\n";
		htmlFile << "     width: 50%; /* Set a width to limit the container size */\n";
		htmlFile << "     margin: 20px auto; /* Center the container and add margin */\n";
		htmlFile << "     padding: 10px; /* Add padding for spacing */\n";
		htmlFile << "     background-color: #f0f0f0; /* Set a background color */\n";
		htmlFile << "     border: 1px solid #ddd; /* Add border for a neat look */\n";
		htmlFile << "     overflow-x: auto; /* Enable horizontal scrolling if needed */\n";
		htmlFile << "     text-align: left;\n";
		htmlFile << "  }\n";
		htmlFile << "  .command-container pre {\n";
		htmlFile << "    white-space: pre-wrap; /* Allow the text to wrap within <pre> */\n";
		htmlFile << "  }\n";
		htmlFile << "    p {\n";
		htmlFile << "      text-align: left; \n";
//		htmlFile << "      max-width: 1100px; \n";
		htmlFile << "      margin: 0, auto; \n";
		htmlFile << "	   padding: 10px; \n";
		htmlFile << "    }\n";
		htmlFile << "    h1 {\n";
		htmlFile << "      font-family: 'Arial', sans-serif; /* Set title font family */\n";
		htmlFile << "      font-size: 36px; /* Set title font size */\n";
		htmlFile << "      margin: 20px 0; /* Add some vertical space */\n";
		htmlFile << "      margin-left: 20px; /* Add some vertical space */\n";
		htmlFile << "      text-align: center; /* Center align title */\n";
		htmlFile << "      text-transform: none; /* Convert text to uppercase */\n";
		htmlFile << "    }\n";
		htmlFile << "    h2 {\n";
		htmlFile << "      font-family: 'Arial', sans-serif; /* Set title font family */\n";
		htmlFile << "      font-size: 24px; /* Set title font size */\n";
		htmlFile << "      color: #555; /* Set title color */\n";
		htmlFile << "      margin: 10px 0; /* Add some vertical space */\n";
		htmlFile << "      margin-left: 2px; /* Add some vertical space */\n";
		htmlFile << "      text-align: left; /* Center align title */\n";
		htmlFile << "      text-transform: none; /* Convert text to uppercase */\n";
		htmlFile << "    }\n";
		htmlFile << "    h3 {\n";
		htmlFile << "      font-family: 'Arial', sans-serif; /* Set title font family */\n";
//		htmlFile << "      font-size: 20px; /* Set title font size */\n";
		htmlFile << "      color: #555; /* Set title color */\n";
		htmlFile << "      margin: 5px 0; /* Add some vertical space */\n";
		htmlFile << "      margin-left: 30px; /* Add some vertical space */\n";
		htmlFile << "      text-align: center; /* Center align title */\n";
		htmlFile << "      text-transform: none; /* Convert text to uppercase */\n";
		htmlFile << "    }\n";
		// Add CSS styles for tables
		htmlFile << "    table {\n";
		htmlFile << "      width: 50%; /* Set the table width to 100% */\n";
		htmlFile << "      border-collapse: collapse;\n";
		htmlFile << "      margin: 20px auto; /* Set table margins */\n"; //auto
		htmlFile << "      margin-bottom: 20px; /* Add some space between tables */\n";
		htmlFile << "      background-color: #f5f5f5; /* Set table background color */\n";
		htmlFile << "    }\n";
		htmlFile << "    th, td {\n";
		htmlFile << "      border: 1px solid #ddd; /* Set border color */\n";
		htmlFile << "      padding: 12px; /* Set cell padding */\n";
		htmlFile << "      text-align: center;\n";
		htmlFile << "    }\n";
		htmlFile << "    th {\n";
		htmlFile << "      background-color: #4CAF50; /* Set header background color */\n";
		htmlFile << "      color: white; /* Set header text color */\n";
		htmlFile << "    }\n";
		htmlFile << "  </style>\n";
		htmlFile << "</head>\n";
		htmlFile << "<body>\n";

		htmlFile << "<h1>Evaluation Reports</h1>\n";
		htmlFile << "<div class=\"command-container\">\n";
		htmlFile << "<pre>SV_STAT command: " << "<br>$ " << program_cmd_str << "</pre>\n";
		htmlFile << "</div>\n";

		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>1. Evaluation results</h2>";
		if(typeMatchLevel.compare(MATCHLEVEL_L)==0){
			htmlFile <<"</div>\n";
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h3 style=\"margin-top: 0px; margin-bottom: 0px; text-align: left;\">Variant type match mode: <span style=\"color:red;\">" << typeMatchLevel<< " (allow type match between DUPLICTION and INSERTION)</span></h3>";
			htmlFile << "</div>";
		}else{
			htmlFile <<"</div>\n";
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h3 style=\"margin-top: 0px; margin-bottom: 0px; text-align: left;\">Variant type match mode: <span style=\"color:red;\">" << typeMatchLevel<< " (strict variant type match)</span></h3>";
			htmlFile << "</div>";
		}
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "The evaluation metrics has two categories after filtering long SV regions: one category is used to highlight performance by metrics including Recall, Precision, F1 score, and sequence consistency (Seqcons) and the other category presents benchmark results, which consists of TP_bench, TP_user, FP, FN. Visualizing these metrics through bar charts provides a more intuitive representation of the assessment results for the variation detection methods.";
		htmlFile << "</div>";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 0px; text-align: left;\">(1) The evaluation results of the user-called set are as follows:</h3>";
		htmlFile << "</div>\n";
		// Add h2 style for the first title
		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Table "<< table_num  <<" Structural Variation Detection Method Performance Evaluation</h4>\n";
		htmlFile << "</div>\n";
		++table_num;
		// Add first table
		htmlFile << "<table border=\"1\">\n";
		// Add header line
		htmlFile << "  <tr>\n";
		for (const auto& title : allmetrics[0]) {
			htmlFile << "    <th>" << title << "</th>\n";
		}
		htmlFile << "  </tr>\n";
		// Traverse the two-dimensional vector a, skipping the first row (header row)
		for (size_t i = 1; i < allmetrics.size(); ++i) {
			htmlFile << "  <tr>\n";
			for (const auto& element : allmetrics[i]) {
				htmlFile << "    <td>" << element << "</td>\n";
			}
			htmlFile << "  </tr>\n";
		}
		htmlFile << "</table>\n";

		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "The table 1 shows the evaluation results of the variation identification result. Where #SVs_bench represents the number of identified structural variations (SVs) in the benchmark set, #SV_user represents the number of SVs in the called set, and #SV_filtered_user represents the number of SVs after filtering out large SVs. #TP stands for the number of True Positives, indicating correctly identified targets or events. #FP stands for the number of False Positives, representing falsely identified targets or events. #FN represents the number of False Negatives, referring to the targets or events that were missed or not identified correctly. Seqcons represents the sequence consistency, which refers to calculating the sequence consistency score for matched SV pairs that include sequences.\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 0px; text-align: left;\">(2) The evaluation results of two categorizes of metrics are shown in the figure:</h3>\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "Two categories of metrics are independently calculated: (a) one category includes Recall, Precision, F1 Score, and Seqcons; (b) the other category consists of #TP_bench, #TP_user, #FP, and #FN. The result statistics are as follows:\n";
		htmlFile << "</div>\n";
//		htmlFile << "<h3>Figure 1. Evaluation results of the user-call set</h3>\n";
		htmlFile << "<style>\n";
		htmlFile << "    .image-container {\n";
		htmlFile << "    display: flex; /* Use flexbox for layout */\n";
		htmlFile << "    justify-content: center; /* Center align images */\n";
		htmlFile << "    margin: 20px auto; /* Add some margin for spacing */\n";
		htmlFile << "  }\n";
		htmlFile << "  .image-container img {\n";
		htmlFile << "    width: 480px; /* Set width of each image */\n";
		htmlFile << "    height: 384px; /* Set height of each image */\n";
		htmlFile << "    border: 1px solid #ddd; /* Add border for a neat look */\n";
		htmlFile << "    border-radius: 5px; /* Add border radius for a rounded look */\n";
		htmlFile << "    margin: 0 10px; /* Add horizontal margin between images */\n";
		htmlFile << "  }\n";
		htmlFile << "</style>\n";
		// Add image container
		htmlFile << "<div class=\"image-container\">\n";
		htmlFile <<"<h4 style=\" margin-left: 40px;\">(a)</h4>\n";
		htmlFile << "  <img src=\""<< figuresFilePath  << "/evaluation_result.png\" alt=\"Benchmark results between different tools\">\n";
		htmlFile <<"<h4 style=\" margin-left: 10px;\">(b)</h4>\n";
		htmlFile << "  <img src=\""<< figuresFilePath  << "/result_classification.png\" alt=\"Performance comparison between different tools\">\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" Evaluation results of the user-call set</h4>\n";
		htmlFile <<"</div>\n";
		++figure_num;

		// Add h2 style for the second title
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>2. Statistical results of deviations for overlapping variants</h2>\n";
		htmlFile <<"</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "For variations that overlap between the user-called set and the benchmark set, the deviations between them are quantified by calculating the center distance and the region size ratio of the overlapping variations.<br>";
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">(1) Deviation of the center distance</h3>";
		htmlFile << "As the center distance approaches 0, the deviation decreases, indicating a more precise identification result. Statistics results for eight size regions are presented in Table 2:";
		htmlFile << "</div>\n";
		htmlFile << "<h3>Table "<< table_num  <<" Statistical results of center distance deviation</h3>\n";
		++table_num;
		// Add second table
		htmlFile << "<table border=\"1\">\n";
		// Add header line
		htmlFile << "  <tr>\n";
		for (const auto& title : centerDistances[0]) {
			htmlFile << "    <th>" << title << "</th>\n";
		}
		htmlFile << "  </tr>\n";
		// Traverse the two-dimensional vector b, skipping the first row (header row)
		for (size_t i = 1; i < centerDistances.size(); ++i) {
			htmlFile << "  <tr>\n";
			for (const auto& element : centerDistances[i]) {
				htmlFile << "    <td>" << element << "</td>\n";
			}
			htmlFile << "  </tr>\n";
		}
		htmlFile << "</table>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">(2) Deviation of the region size ratio</h3>";
		htmlFile << "Calculating the region size ratio for two overlapping variations based on the length of SVs, the closer the ratio is to 1, the smaller the deviation, indicating a more precise and accurate identification result. Statistics results for nine size regions are presented in Table 3:";
		htmlFile << "</div>\n";
		// Add h2 style for the third title
		htmlFile << "<h3>Table "<< table_num  <<" Statistical results of deviation of the region size ratio</h3>\n";
		++table_num;
		// Add third table
		htmlFile << "<table border=\"1\">\n";
		// Add header line
		htmlFile << "  <tr>\n";
		for (const auto& title : sizeratios[0]) {
			htmlFile << "    <th>" << title << "</th>\n";
		}
		htmlFile << "  </tr>\n";
		// Traverse the two-dimensional vector c, skipping the first row (header row)
		for (size_t i = 1; i < sizeratios.size(); ++i) {
			htmlFile << "  <tr>\n";
			for (const auto& element : sizeratios[i]) {
				htmlFile << "    <td>" << element << "</td>\n";
			}
			htmlFile << "  </tr>\n";
		}
		htmlFile << "</table>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 0px;\">The statistical results of user-called sets are as follows:</h3>";
		htmlFile << "</div>\n";
		if(folderPng3.size()==1){
//			htmlFile << "<div style=\"text-align: center; margin: 0 auto;\">\n";
			htmlFile << "<div style=\"display: flex; flex-wrap: wrap; justify-content: center; margin: 0 auto;\">\n";
			htmlFile << "  <img src=\""<< folderPng3[0]  << "\" alt=\"result details\" style=\"width: 40%;\">\n";
			htmlFile << "</div>\n";
		}else{
			/*htmlFile << "<div style=\"display: flex; flex-wrap: wrap; justify-content: center; margin: 0 auto;\">\n";
			for (size_t i = 0; i < folderPng3.size(); ++i) {
//				htmlFile << "  <div style=\"text-align: center;\">\n";
				htmlFile <<"<h4 style=\"margin-top: 5px; margin-bottom: 5px; display:block;\">(" << i+1 << ") " << allmetrics[i+1][0] << ":</h4>\n";
				htmlFile << "  <img src=\"" << folderPng3[i] << "\" alt=\"result details" << i + 1 << "\"style=\"width: 32%; margin: 5px;\">\n";
//				htmlFile << "</div>\n";
				// Add a line break after every two images
				if ((i + 1) % 2 == 0) {
					htmlFile << "<br>\n";
				}
			}
			htmlFile << "</div>\n";*/
			htmlFile << "<div style=\"display: flex; flex-wrap: wrap; justify-content: center; margin: 0 auto;\">\n";
//			size_t asciiValue = static_cast<int>('a');
			for (size_t i = 0, asciiValue = static_cast<int>('a'); i < folderPng3.size(); i += 2) {
				htmlFile << "<div style=\"display: flex; justify-content: center;\">\n";
				if(i == folderPng3.size()-1 and folderPng3.size()%2 ==1){
					htmlFile << "<div style=\"display: flex; flex-wrap: wrap; justify-content: center; margin: 0 auto;\">\n";
					htmlFile << "<h4 style=\"text-align:center; margin-top: 0px; margin-bottom: 5px;\">(" << static_cast<char>(asciiValue+ i) << ") " << allmetrics[i+1][0] << ":</h4>\n";
					htmlFile << "<img src=\""<< folderPng3[i]  << "\" alt=\"result details\" style=\"width: 45%; margin: 10px;\">\n";
					htmlFile << "</div>\n";
					htmlFile << "</div>\n";
					break;
				}
				// Image and text for the first image in the row
				htmlFile << "<div style=\"text-align: right; margin-right: 2px;\">\n";
				htmlFile << "<h4 style=\"text-align:center; margin-top: 0px; margin-bottom: 5px;\">(" << static_cast<char>(asciiValue+ i) << ") " << allmetrics[i+1][0] << ":</h4>\n";
				htmlFile << "<img src=\"" << folderPng3[i] << "\" alt=\"result details " << i + 1 << "\" style=\"width: 70%; margin: 5px;\">\n";
				htmlFile << "</div>\n";

				// Check if there is a second image in the row
				if (i + 1 < folderPng3.size()) {
					// Image and text for the second image in the row
					htmlFile << "<div style=\"text-align: left; margin-left: 2px;\">\n";
					htmlFile << "<h4 style=\"text-align:left; margin-top: 0px; margin-bottom: 5px;\">(" << static_cast<char>(asciiValue + i + 1) << ") " << allmetrics[i+2][0] << ":</h4>\n";
					htmlFile << "<img src=\"" << folderPng3[i + 1] << "\" alt=\"result details " << i + 2 << "\" style=\"width: 70%; margin: 5px;\">\n";
					htmlFile << "</div>\n";
				}

				htmlFile << "</div>\n";
			}
			htmlFile << "</div>\n";
		}
		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" Deviation statistics with overlapping variations</h4>\n";
		htmlFile <<"</div>\n";
		++figure_num;

		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>3. Evaluation results for metrics of different SV size regions</h2>\n";
		htmlFile <<"</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "The SV identification results typically contain variations of various sizes, and categorize these variations into different size ranges could be used to explore the identification results more detailed in a fine-grained manner, and could provide new insights into the sensitivity of SV callers to variations of different sizes. Detailed evaluation results are presented in the table as follows：\n";
		htmlFile << "</div>\n";
		if(folderPng5.size()>0){
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">(1) Evaluation results for metrics of different SV size regions with different methods</h3>";
			htmlFile <<"Variations are categorized into seven size regions and metrics are computed for comprehensive evaluation for different detection methods within each region. The evaluation results are as follows:\n";
			htmlFile << "</div>\n";
			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
				htmlFile << "  <img src=\""<< folderPng5[0] << "\" alt=\"different_range.png\" style=\"width: 50%;\">\n";
			htmlFile << "</div>\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" Statistics of metrics of different SV size region</h4>\n";
			htmlFile <<"</div>\n";
			++figure_num;
		}
		//The metric evaluation results of cutesv in different SV regions
		for (size_t tableIndex = 0; tableIndex < allregionmetrics.size(); ++tableIndex) {
//			("<<tableIndex+1 <<")
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			if(folderPng5.size()>0){
				htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">("<<tableIndex+ 2 <<") The user-called set ("<< alltoolnames[tableIndex] <<") of basic metrics results statistics</h3>";
			}else{
				htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">("<<tableIndex+1 <<") The user-called set ("<< alltoolnames[tableIndex] <<") of basic metrics results statistics</h3>";
			}
			htmlFile << "</div>\n";
			// Add h2 style for the title of each table
			htmlFile << "<h3> Table "<< table_num  <<"  The metric evaluation results of " << alltoolnames[tableIndex] << " in different SV regions</h3>\n";
			++table_num;
			htmlFile << "<table border=\"1\">\n";
			// Add header line
			htmlFile << "  <tr>\n";
			for (const auto& title : allregionmetrics[tableIndex][0]) {
				htmlFile << "    <th>" << title << "</th>\n";
			}
			htmlFile << "  </tr>\n";

			// Traverse the two-dimensional vector d[tableIndex], skipping the first row (header row)
			for (size_t i = 1; i < allregionmetrics[tableIndex].size(); ++i) {
				htmlFile << "  <tr>\n";
				for (const auto& element : allregionmetrics[tableIndex][i]) {
					htmlFile << "    <td>" << element << "</td>\n";
				}
				htmlFile << "  </tr>\n";
			}
			htmlFile << "</table>\n";
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4>Evaluation results for metrics of different SV size regions show as following figures:</h4>";
			htmlFile << "</div>\n";
			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
			htmlFile <<"<h4 style=\" margin-left: 80px;\">(a)</h4>\n";
			htmlFile << "  <img src=\""<< folderPng4[2*tableIndex] << "\" alt=\" calculation_result.png\" style=\"width: 28%;\">\n";
			htmlFile <<"<h4 style=\" margin-left: 5px;\">(b)</h4>\n";
			htmlFile << "  <img src=\""<< folderPng4[2*tableIndex+1]  << "\" alt=\" quantitative_Statistics.png\" style=\"width: 28%;\">\n";
			htmlFile << "</div>\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" Result statistics of different SV size region</h4>\n";
			htmlFile <<"</div>\n";
			++figure_num;
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "Figure (a) shows the statistical results of Recall, Precision, F1 score and Seqcons; (b) shows the statistical results of #TP_benchmark, #TP_user,#FP and #FN.";
			htmlFile <<"</div>\n";
		}

		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>4. Quantitative statistics</h2>\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		//htmlFile << "The benchmark set：<br>";
		htmlFile << "<h4 style=\"margin-top: 0px; margin-bottom: 5px;\">(1) Statistics of the count of different SV lengths in the benchmark set: </h4>";
		htmlFile << "The SV reference region size statistics for benchmark set: Total SVs number："<< allmetrics[1][1] << "";
		htmlFile << "<div style=\"text-align: center; margin: 0 auto;\">\n";
		htmlFile << "  <img src=\""<< folderPng1[0]  << "\" alt=\"The benchmark set SVs count statistics\" style=\"width: 60%;\">\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" The quantity distribution of the benchmark set</h4>\n";
		htmlFile <<"</div>\n";
		++figure_num;
		htmlFile <<"The figure shows the distribution of SV counts of the benchmark set.\n";
		htmlFile <<"<br>\n";
		//htmlFile << "The user-called set：<br>";
		for(size_t i=1; i<allmetrics.size(); i++){
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 5px;\">("<< i+1 <<") Statistics of the count of different SV lengths in the user-called set ("<< allmetrics[i][0] <<"): " << "</h4>";
			htmlFile << "The SV reference region size statistics before filtering for user-called set ("<< allmetrics[i][0]<<"):"<< "Total SVs number：" << allmetrics[i][2]<<"<br>";
			htmlFile << "The SV reference region size statistics after filtering for user-called set ("<< allmetrics[i][0]<<"):"<< "Total SVs number：" << allmetrics[i][3]<<"<br>";
			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
			htmlFile <<"<h4 style=\" margin-left: 40px;\">(a)</h4>\n";
			htmlFile << "  <img src=\""<< folderPng1[2*(i-1)+1] << "\" alt=\" The user-called set SVs count statistics\" style=\"width: 50%;\">\n";
			htmlFile <<"<h4 style=\" margin-left: 10px;\">(b)</h4>\n";
			htmlFile << "  <img src=\""<< folderPng1[2*(i-1)+2]  << "\" alt=\" The user-called set filtered SVs count statistics\" style=\"width: 50%;\">\n";
			htmlFile << "</div>\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" The quantity distribution of the user-called set</h4>\n";
			htmlFile <<"</div>\n";
			++figure_num;
			htmlFile << "The figures show the distribution of SV counts, where (a) represents the result statistics before filtering large SVs, and (b) shows the result statistics after filtering large SVs.\n";
		}
		htmlFile << "</div>\n";

		htmlFile << "<br>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>More information</h2>\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<ul>\n";
		htmlFile << "<li>For more detailed evaluation results, please refer to the generated result information in the respective folders.</li>\n";
		htmlFile << "<li>For more detailed experiment information, please refer to the github repositories: " << " <a href=\"https://github.com/zhuxiao/sv_stat\" target=\"_blank\">sv_stat</a>" <<  " and <a href=\"https://github.com/zhuxiao/sv_stat-experiments\" target=\"_blank\">sv_stat-experiments</a>.</li>\n";
		htmlFile << "<li>If you have any problems, comments, or suggestions, please contact xzhu@ytu.edu.cn without hesitation. Thank you very much!</li>\n";
		htmlFile << "</ul>\n";
		htmlFile << "-------------------------- This is the end of the Evaluation Reports. --------------------------\n";
		htmlFile << "</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h2>2. Evaluation results</h2>";
//		htmlFile <<"</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "Categorizing the evaluation metrics into two categories after filtering long SV regions: one highlighting performance metrics including recall, precision, and F1 score, and the other presenting benchmark results, encompassing TP_bench, TP_user, FP, FN.   Visualizing these metrics through bar charts offers a more intuitive representation of the assessment results for the detection method.";
//		htmlFile << "</div>";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">The evaluation results of the user-called set are as follows:</h4>";
//		htmlFile << "</div>\n";
//		// Add h2 style for the first title
//		htmlFile << "<h3>Table 1 Structural Variation Detection Method Performance Evaluation</h3>\n";
//		// Add first table
//		htmlFile << "<table border=\"1\">\n";
//		// Add header line
//		htmlFile << "  <tr>\n";
//		for (const auto& title : allmetrics[0]) {
//			htmlFile << "    <th>" << title << "</th>\n";
//		}
//		htmlFile << "  </tr>\n";
//		// Traverse the two-dimensional vector a, skipping the first row (header row)
//		for (size_t i = 1; i < allmetrics.size(); ++i) {
//			htmlFile << "  <tr>\n";
//			for (const auto& element : allmetrics[i]) {
//				htmlFile << "    <td>" << element << "</td>\n";
//			}
//			htmlFile << "  </tr>\n";
//		}
//		htmlFile << "</table>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "The table 1 shows the evaluation results of the structural variation detection method. Where SVs_bench represents the quantity of structural variations (SVs) in the benchmark set, SV_user represents the quantity of SVs in the called set, and SV_user_filtered represents the quantity of SVs after filtering out large SVs. TP stands for True Positive, indicating the count of correctly identified targets or events. FP stands for False Positive, representing the count of falsely identified targets or events. FN represents False Negative, referring to the count of targets or events that were missed or not identified correctly.\n";
//		htmlFile << "</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">The evaluation results of two categorizes of metrics are shown in the figure:</h4>\n";
//		htmlFile << "</div>\n";
////		htmlFile << "<h3>Figure 1. Evaluation results of the user-call set</h3>\n";
//		htmlFile << "<style>\n";
//		htmlFile << "    .image-container {\n";
//		htmlFile << "    display: flex; /* Use flexbox for layout */\n";
//		htmlFile << "    justify-content: center; /* Center align images */\n";
//		htmlFile << "    margin: 20px auto; /* Add some margin for spacing */\n";
//		htmlFile << "  }\n";
//		htmlFile << "  .image-container img {\n";
//		htmlFile << "    width: 500px; /* Set width of each image */\n";
//		htmlFile << "    height: 400px; /* Set height of each image */\n";
//		htmlFile << "    border: 1px solid #ddd; /* Add border for a neat look */\n";
//		htmlFile << "    border-radius: 5px; /* Add border radius for a rounded look */\n";
//		htmlFile << "    margin: 0 10px; /* Add horizontal margin between images */\n";
//		htmlFile << "  }\n";
//		htmlFile << "</style>\n";
//		// Add image container
//		htmlFile << "<div class=\"image-container\">\n";
//		htmlFile << "  <img src=\""<< figuresFilePath  << "/evaluation_result.png\" alt=\"Benchmark results between different tools\">\n";
//		htmlFile << "  <img src=\""<< figuresFilePath  << "/result_classification.png\" alt=\"Performance comparsion between different tools\">\n";
//		htmlFile << "</div>\n";
//		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure 1. Evaluation results of the user-call set</h4>\n";
//		htmlFile <<"</div>\n";

//		// Add h2 style for the second title
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h2>3. Statistical results of deviations for overlapping variants</h2>\n";
//		htmlFile <<"</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "For variations that overlap between the user-called set and the benchmark set, the deviations between them are quantified by calculating the center distance and the region size ratio of the overlapping variations.<br>";
//		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">(1) Deviation of the center distance</h3>";
//		htmlFile << "As the center distance approaches 0, the deviation decreases, indicating a more precise identification result. Statistics results for eight size regions are presented in Table 2:";
//		htmlFile << "</div>\n";
//		htmlFile << "<h3>Table 1 Statistical results of center distance deviation</h3>\n";
//		// Add second table
//		htmlFile << "<table border=\"1\">\n";
//		// Add header line
//		htmlFile << "  <tr>\n";
//		for (const auto& title : centerDistances[0]) {
//			htmlFile << "    <th>" << title << "</th>\n";
//		}
//		htmlFile << "  </tr>\n";
//		// Traverse the two-dimensional vector b, skipping the first row (header row)
//		for (size_t i = 1; i < centerDistances.size(); ++i) {
//			htmlFile << "  <tr>\n";
//			for (const auto& element : centerDistances[i]) {
//				htmlFile << "    <td>" << element << "</td>\n";
//			}
//			htmlFile << "  </tr>\n";
//		}
//		htmlFile << "</table>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">(2) Deviation of the region size ratio</h3>";
//		htmlFile << "Calculating the region size ratio for two overlapping variations based on the length of SVs, the closer the ratio is to 1, the smaller the deviation, indicating a more precise and accurate identification result. Statistics results for nine size regions are presented in Table 3:";
//		htmlFile << "</div>\n";
//		// Add h2 style for the third title
//		htmlFile << "<h3>Table 2 Statistical results of deviation of the region size ratio</h3>\n";
//		// Add third table
//		htmlFile << "<table border=\"1\">\n";
//		// Add header line
//		htmlFile << "  <tr>\n";
//		for (const auto& title : sizeratios[0]) {
//			htmlFile << "    <th>" << title << "</th>\n";
//		}
//		htmlFile << "  </tr>\n";
//		// Traverse the two-dimensional vector c, skipping the first row (header row)
//		for (size_t i = 1; i < sizeratios.size(); ++i) {
//			htmlFile << "  <tr>\n";
//			for (const auto& element : sizeratios[i]) {
//				htmlFile << "    <td>" << element << "</td>\n";
//			}
//			htmlFile << "  </tr>\n";
//		}
//		htmlFile << "</table>\n";
//		htmlFile << "<div style=\"text-align: center; margin: 0 auto;\">\n";
//		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 0px;\">The statistical results of user-called sets are as followed:</h3>";
//		htmlFile << "</div>\n";
//		if(folderPng3.size()==1){
////			htmlFile << "<div style=\"text-align: center; margin: 0 auto;\">\n";
//			htmlFile << "<div style=\"display: flex; flex-wrap: wrap; justify-content: center; margin: 0 auto;\">\n";
//			htmlFile << "  <img src=\""<< folderPng3[0]  << "\" alt=\"result details\" style=\"width: 50%;\">\n";
//			htmlFile << "</div>\n";
//		}else{
//			/*htmlFile << "<div style=\"display: flex; flex-wrap: wrap; justify-content: center; margin: 0 auto;\">\n";
//			for (size_t i = 0; i < folderPng3.size(); ++i) {
////				htmlFile << "  <div style=\"text-align: center;\">\n";
//				htmlFile <<"<h4 style=\"margin-top: 5px; margin-bottom: 5px; display:block;\">(" << i+1 << ") " << allmetrics[i+1][0] << ":</h4>\n";
//				htmlFile << "  <img src=\"" << folderPng3[i] << "\" alt=\"result details" << i + 1 << "\"style=\"width: 32%; margin: 5px;\">\n";
////				htmlFile << "</div>\n";
//				// Add a line break after every two images
//				if ((i + 1) % 2 == 0) {
//					htmlFile << "<br>\n";
//				}
//			}
//			htmlFile << "</div>\n";*/
//			htmlFile << "<div style=\"display: flex; flex-wrap: wrap; justify-content: center; margin: 0 auto;\">\n";
////			size_t asciiValue = static_cast<int>('a');
//			for (size_t i = 0, asciiValue = static_cast<int>('a'); i < folderPng3.size(); i += 2) {
//			    htmlFile << "<div style=\"display: flex; justify-content: center;\">\n";
//			    if(i == folderPng3.size()-1 and folderPng3.size()%2 ==1){
//			    	htmlFile << "<div style=\"display: flex; flex-wrap: wrap; justify-content: center; margin: 0 auto;\">\n";
//			    	htmlFile << "<h4 style=\"text-align:center; margin-top: 0px; margin-bottom: 5px;\">(" << static_cast<char>(asciiValue+ i) << ") " << allmetrics[i+1][0] << ":</h4>\n";
//					htmlFile << "<img src=\""<< folderPng3[i]  << "\" alt=\"result details\" style=\"width: 40%; margin: 10px;\">\n";
//					htmlFile << "</div>\n";
//					htmlFile << "</div>\n";
//					break;
//			    }
//			    // Image and text for the first image in the row
//			    htmlFile << "<div style=\"text-align: right; margin-right: 2px;\">\n";
//			    htmlFile << "<h4 style=\"text-align:center; margin-top: 0px; margin-bottom: 5px;\">(" << static_cast<char>(asciiValue+ i) << ") " << allmetrics[i+1][0] << ":</h4>\n";
//			    htmlFile << "<img src=\"" << folderPng3[i] << "\" alt=\"result details " << i + 1 << "\" style=\"width: 70%; margin: 5px;\">\n";
//			    htmlFile << "</div>\n";
//
//			    // Check if there is a second image in the row
//			    if (i + 1 < folderPng3.size()) {
//			        // Image and text for the second image in the row
//			        htmlFile << "<div style=\"text-align: left; margin-left: 2px;\">\n";
//			        htmlFile << "<h4 style=\"text-align:left; margin-top: 0px; margin-bottom: 5px;\">(" << static_cast<char>(asciiValue + i + 1) << ") " << allmetrics[i+2][0] << ":</h4>\n";
//			        htmlFile << "<img src=\"" << folderPng3[i + 1] << "\" alt=\"result details " << i + 2 << "\" style=\"width: 70%; margin: 5px;\">\n";
//			        htmlFile << "</div>\n";
//			    }
//
//			    htmlFile << "</div>\n";
//			}
//			htmlFile << "</div>\n";
//		}
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h2>4. Evaluation results for metrics of different SV size regions</h2>\n";
//		htmlFile <<"</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "The SV identification results typically contain variations of various sizes, and categorize these variations into different size ranges could be used to explore the identification results more detailed in a fine-grained manner, and could provide new insights into the sensitivity of SV callers to variations of different sizes. Detailed evaluation results are presented in the table as follows：\n";
//		htmlFile << "</div>\n";
//		//The metric evaluation results of cutesv in different SV regions
//		for (size_t tableIndex = 0; tableIndex < allregionmetrics.size(); ++tableIndex) {
////			("<<tableIndex+1 <<")
//			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//			htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">("<<tableIndex+1 <<") The user-called set ("<< alltoolnames[tableIndex] <<") of basic metrics results statistics</h3>";
//			htmlFile << "</div>\n";
//			// Add h2 style for the title of each table
//			htmlFile << "<h3> Table "<<(size_t)1 + tableIndex << " The metric evaluation results of " << alltoolnames[tableIndex] << " in different SV regions</h3>\n";
//			htmlFile << "<table border=\"1\">\n";
//			// Add header line
//			htmlFile << "  <tr>\n";
//			for (const auto& title : allregionmetrics[tableIndex][0]) {
//				htmlFile << "    <th>" << title << "</th>\n";
//			}
//			htmlFile << "  </tr>\n";
//
//			// Traverse the two-dimensional vector d[tableIndex], skipping the first row (header row)
//			for (size_t i = 1; i < allregionmetrics[tableIndex].size(); ++i) {
//				htmlFile << "  <tr>\n";
//				for (const auto& element : allregionmetrics[tableIndex][i]) {
//					htmlFile << "    <td>" << element << "</td>\n";
//				}
//				htmlFile << "  </tr>\n";
//			}
//			htmlFile << "</table>\n";
//			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//			htmlFile << "<h4>Evaluation results for metrics of different SV size regions show as follow with figures:</h4>";
//			htmlFile << "</div>\n";
//			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
//				htmlFile << "  <img src=\""<< folderPng4[2*tableIndex] << "\" alt=\" calculation_result.png\" style=\"width: 30%;\">\n";
//				htmlFile << "  <img src=\""<< folderPng4[2*tableIndex+1]  << "\" alt=\" quantitative_Statistics.png\" style=\"width: 30%;\">\n";
//			htmlFile << "</div>\n";
//		}
//
//		if(folderPng5.size()>0){
//			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//			htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">("<< allregionmetrics.size() +1 <<") Evaluation results for metrics of different SV size regions with different methods</h3>";
//			htmlFile <<"Variations are categorized into seven size regions and metrics are computed for comprehensive evaluation for different detection methods within each region. The evaluation results are as follows:\n";
//			htmlFile << "</div>\n";
//			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
//				htmlFile << "  <img src=\""<< folderPng5[0] << "\" alt=\"different_range.png\" style=\"width: 50%;\">\n";
//			htmlFile << "</div>\n";
//		}
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "For more detailed evaluation results, please refer to the generated result information in the respective folders.<br>\n";
//		htmlFile << "</div>\n";
//		htmlFile << "<h2>5. Evaluation Results with Bar Chart </h2>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile <<"Categorizing the evaluation metrics into two categories: one highlighting performance metrics including recall, precision, and F1 score, and the other presenting benchmark results, encompassing TP_bench, TP_user, FP, FN.   Visualizing these metrics through bar charts offers a more intuitive representation of the assessment results for the detection method. The evaluation results are shown in the figure:\n";
//		htmlFile << "</div>\n";
//		htmlFile << "<h3>Comprehensive evaluation results</h3>\n";
//		htmlFile << "<style>\n";
//		htmlFile << "    .image-container {\n";
//		htmlFile << "    display: flex; /* Use flexbox for layout */\n";
//		htmlFile << "    justify-content: center; /* Center align images */\n";
//		htmlFile << "    margin: 20px auto; /* Add some margin for spacing */\n";
//		htmlFile << "  }\n";
//		htmlFile << "  .image-container img {\n";
//		htmlFile << "    width: 500px; /* Set width of each image */\n";
//		htmlFile << "    height: 400px; /* Set height of each image */\n";
//		htmlFile << "    border: 1px solid #ddd; /* Add border for a neat look */\n";
//		htmlFile << "    border-radius: 5px; /* Add border radius for a rounded look */\n";
//		htmlFile << "    margin: 0 10px; /* Add horizontal margin between images */\n";
//		htmlFile << "  }\n";
//		htmlFile << "</style>\n";
//		// Add image container
//		htmlFile << "<div class=\"image-container\">\n";
//		htmlFile << "  <img src=\""<< figuresFilePath  << "/evaluation_result.png\" alt=\"Benchmark results between different tools\">\n";
//		htmlFile << "  <img src=\""<< figuresFilePath  << "/result_classification.png\" alt=\"Performance comparsion between different tools\">\n";
//		htmlFile << "</div>\n";

		// End HTML file
	   htmlFile << "</body>\n";
	   htmlFile << "</html>\n";
	   htmlFile.close();

	  }else {
	   cerr << "Unable to open HTML file.\n";
	   }
}
