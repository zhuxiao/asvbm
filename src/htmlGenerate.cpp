#include "htmlGenerate.h"


std::string escapeHtml(const std::string& str) {
    std::ostringstream escaped;
    for (char ch : str) {
        switch (ch) {
            case '&':  escaped << "&amp;";    break;
            case '\"': escaped << "&quot;";   break;
            case '\'': escaped << "&#39;";    break;
            case '<':  escaped << "&lt;";     break;
            case '>':  escaped << "&gt;";     break;
            default:   escaped << ch;        break;
        }
    }
    return escaped.str();
}

void Generatehtml(string figuresFilePath){
	int32_t table_num = 1, figure_num = 1;
	string HtmlFilePath = outputPathname + '/' + htmlFilename;
	ofstream htmlFile(HtmlFilePath);

//	std::string escaped_cmd_str = escapeHtml(all_cmd_program);

	if (htmlFile.is_open()) {
		htmlFile << "<!DOCTYPE html>\n";
		htmlFile << "<html>\n";
		htmlFile << "<head>\n";
		htmlFile << "  <title>Benchmarking Reports</title>\n";
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

		htmlFile << "<h1>Benchmarking Reports</h1>\n";
		htmlFile << "<div class=\"command-container\"  style=\"max-width: 100%; overflow-x: auto;\">\n";
		htmlFile << "<pre style=\"overflow-x: auto; white-space: pre-wrap; word-wrap: break-word;\">ASVBM command: " << "<br>$ " << all_cmd_program << "</pre>\n";
		htmlFile << "</div>\n";

		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>1. Benchmarking results</h2>";
		if(typeMatchLevel.compare(MATCHLEVEL_L)==0){
			htmlFile <<"</div>\n";
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h3 style=\"margin-top: 0px; margin-bottom: 0px; text-align: left;\">ASVBM version: <span style=\"color:red;\">" << PROG_VERSION<< "</span></h3>";
			htmlFile << "<h3 style=\"margin-top: 0px; margin-bottom: 0px; text-align: left;\">Variant type match mode: <span style=\"color:red;\">" << typeMatchLevel<< " (allow type match between DUPLICTION and INSERTION)</span></h3>";
			htmlFile << "</div>";
		}else{
			htmlFile <<"</div>\n";
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h3 style=\"margin-top: 0px; margin-bottom: 0px; text-align: left;\">Variant type match mode: <span style=\"color:red;\">" << typeMatchLevel<< " (strict variant type match)</span></h3>";
			htmlFile << "</div>";
		}
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "The benchmarking metrics has two categories after filtering long Structural Variant(SV) regions: one category is used to highlight performance by metrics including Recall, Precision, F1 score, and sequence similarity (SeqSim) and the other category presents benchmark results, which consists of TP, FP, FN and LP. Visualizing these metrics through bar charts provides a more intuitive representation of the benchmarking results for the variant detection methods.";
		htmlFile << "</div>";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 0px; text-align: left;\">(1) The benchmarking results of the user-called set are as follows:</h3>";
		htmlFile << "</div>\n";
		// Add h2 style for the first title
		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Table "<< table_num  <<" Structural variant Detection Method Performance Benchmarking</h4>\n";
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
		htmlFile << "The table 1 shows the benchmarking results of the variant identification result. Where #SVs_bench represents the count of identified structural variants (SVs) in the benchmark set after filtering out long SVs and short SVs. #SVs_user represents the count of SVs in the called set after filtering out long SVs and short SVs. #TP stands for the count of True Positives, indicating correctly identified targets or events. #FP stands for the count of False Positives, representing falsely identified targets or events. #FN represents the count of False Negatives, referring to the targets or events that were missed or not identified correctly. #LP_user represents variant calls that, referring to adjacent or overlapping variants in the user-called set, collectively resemble true positives in the benchmark set. #LP_bench represents variant calls that, referring to adjacent or overlapping variants in the benchmark set, collectively resemble true positives in the user-called set. SeqSim represents the measure of sequence similarity, which is calculated for matched SV pairs that include sequences.\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 0px; text-align: left;\">(2) The benchmarking results of two categories of metrics are shown in the figure:</h3>\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "Two categories of metrics are independently calculated: (a) one category includes Recall, Precision, F1 Score, and SeqSim; (b) the other category consists of #TP, #FP and #FN. The result statistics are as follows:\n";
		htmlFile << "</div>\n";
//		htmlFile << "<h3>Figure 1. Benchmarking results of the user-call set</h3>\n";
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
		htmlFile << "  <img src=\""<< figuresFilePath  << "/benchmarking_result.png\" alt=\"Benchmark results between different tools\">\n";
		htmlFile <<"<h4 style=\" margin-left: 10px;\">(b)</h4>\n";
		htmlFile << "  <img src=\""<< figuresFilePath  << "/result_classification.png\" alt=\"Performance comparison between different tools\">\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" Benchmarking results of the user-call set</h4>\n";
		htmlFile <<"</div>\n";
		++figure_num;
		if(SVcallernames.size()>1){
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 0px; text-align: left;\">(3) The UpSet plot generated by TP_bench for multiple callsets benchmarking is shown as follows:</h3>\n";
			htmlFile << "</div>\n";

			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "The UpSet plot illustrates the benchmarking of TP_bench variants generated by TP_bench across multiple user callsets. The plot displays the distribution and intersection of high-confidence variants within the benchmark set.\n";
			htmlFile << "</div>\n";

	//		htmlFile << "<div style=\"text-align: center; margin: 0 auto;\">\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto; max-width:800px;\">\n";
	//		htmlFile << "  <img src=\""<< RscriptfilenamePath << "\" alt=\"UpSet plot of TP_bench\" style=\"width: 50%; \">\n"; //style=\"width: 100%;
			htmlFile << "  <img src=\"" << RscriptfilenamePath << "\" alt=\"UpSet plot of TP_bench\" style=\"width: 100%; max-width: 800px; height: auto; object-fit: contain;\">\n";
	//		<img src="path/to/image.png" alt="UpSet plot of TP_bench" style="width: 50%; max-width: 800px; height: auto; object-fit: contain;">
			htmlFile << "</div>\n";

			htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" UpSet plot of TP_bench across multiple callsets</h4>\n";
			htmlFile <<"</div>\n";
		++figure_num;
		}
		// Add h2 style for the second title
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>2. Statistical results of deviations for overlapping variants</h2>\n";
		htmlFile <<"</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "For variants that overlap between the user-called set and the benchmark set, the deviations between them are quantified by calculating the breakpoint distance and the variant size ratio of the overlapping variants.<br>";
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">(1) Deviation of the breakpoint distance</h3>";
		htmlFile << "As the breakpoint distance approaches 0, the deviation decreases, indicating a more precise identification result. Statistics results for eight size regions are presented in Table 2:";
		htmlFile << "</div>\n";
		htmlFile << "<h3>Table "<< table_num  <<" Statistical results of breakpoint distance deviation</h3>\n";
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
		htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">(2) Deviation of the variant size ratio</h3>";
		htmlFile << "Calculating the variant size ratio for two overlapping variants based on the length of SVs, the closer the ratio is to 1, the smaller the deviation, indicating a more precise and accurate identification result. Statistics results for nine size regions are presented in Table 3:";
		htmlFile << "</div>\n";
		// Add h2 style for the third title
		htmlFile << "<h3>Table "<< table_num  <<" Statistical results of deviation of the variant region size ratio</h3>\n";
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
		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" Deviation statistics with overlapping variants</h4>\n";
		htmlFile <<"</div>\n";
		++figure_num;

		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>3. Benchmarking results for metrics of different SV size regions</h2>\n";
		htmlFile <<"</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "The SV identification results typically contain variants of various sizes, and categorize these variants into different size ranges could be used to explore the identification results more detailed in a fine-grained manner, and could provide new insights into the sensitivity of SV callers to variants of different sizes. Detailed benchmarking results are presented in the table as follows：\n";
		htmlFile << "</div>\n";
		if(folderPng5.size()>0){
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">(1) Benchmarking results for metrics of different SV size regions with different methods</h3>";
			htmlFile <<"variants are categorized into eight size regions and metrics are computed for comprehensive benchmarking for different detection methods within each region. The benchmarking results are as follows:\n";
			htmlFile << "</div>\n";
			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
				htmlFile << "  <img src=\""<< folderPng5[0] << "\" alt=\"different_range.png\" style=\"width: 50%;\">\n";
			htmlFile << "</div>\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" Statistics of metrics of different SV size region</h4>\n";
			htmlFile <<"</div>\n";
			++figure_num;
		}
		//The metric benchmarking results of cutesv in different SV regions
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
			htmlFile << "<h3> Table "<< table_num  <<"  The metric benchmarking results of " << alltoolnames[tableIndex] << " in different SV regions</h3>\n";
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
			htmlFile << "<h4>Benchmarking results for metrics of different SV size regions show as following figures:</h4>";
			htmlFile << "</div>\n";
			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
			htmlFile <<"<h4 style=\" margin-left: 80px;\">(a)</h4>\n";
			htmlFile << "  <img src=\""<< folderPng4[2*tableIndex] << "\" alt=\" calculation_result.png\" style=\"width: 28%;\">\n";
			htmlFile <<"<h4 style=\" margin-left: 5px;\">(b)</h4>\n";
			htmlFile << "  <img src=\""<< folderPng4[2*tableIndex+1]  << "\" alt=\" quantitative_Statistics.png\" style=\"width: 28%;\">\n";
			htmlFile << "</div>\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" Result statistics of different SV size ranges</h4>\n";
			htmlFile <<"</div>\n";
			++figure_num;
			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "Figure (a) shows the statistical results of Recall, Precision, F1 score and SeqSim; (b) shows the statistical results of #TP_user, #TP_bench, #FP and #FN.";
			htmlFile <<"</div>\n";
		}

		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>4. SV size distribution statistics</h2>\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
			if(SVcallernames.size()>1){
	//		htmlFile << "<h3 style=\"margin-top: 0px; margin-bottom: 5px;\">(I) Distribution of SV  of multiple user callsets</h3>\n";
			htmlFile << "<h4 style=\"font-size: 1.17em; margin-top: 0px; margin-bottom: 5px;\">(I) Distribution of SVs  of multiple user callsets</h4>\n";
	//		htmlFile << "</div>\n";
	//		htmlFile << "<div style=\"text-align: center; margin: 0 auto;\">\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto;\">\n";
			htmlFile << "  <img src=\""<< SVdistributionPath  << "\" alt=\" SV distribution for mutiple user callsets\" style=\"width: 80%;\">\n";
			htmlFile << "</div>\n";
			htmlFile << "</div>\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" SV distribution for mutiple user callsets</h4>\n";
			htmlFile <<"</div>\n";
			++figure_num;

			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
	//		htmlFile << "<h3 style=\"margin-top: 0px; margin-bottom: 5px;\">(II) Distribution of SV  of benchmark set and user callsets</h3>\n";
			htmlFile << "<h4 style=\"font-size: 1.17em; margin-top: 0px; margin-bottom: 5px;\">(II) Distribution of SVs  of benchmark set and user callsets</h4>\n";
//			htmlFile <<"</div>\n";
		}
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		//htmlFile << "The benchmark set：<br>";
		htmlFile << "<h4 style=\"margin-top: 0px; margin-bottom: 5px;\">(1) Statistics of the count of different SV lengths in the benchmark set: </h4>";
		htmlFile << "The SV reference region size statistics before filtering for benchmark set: Total SV count："<< before_filtered[0]<<"<br>";
		htmlFile << "The SV reference region size statistics after filtering for benchmark set: Total SV count："<< allmetrics[1][1] << "";
/*		htmlFile << "<div style=\"text-align: center; margin: 0 auto;\">\n";
		htmlFile << "  <img src=\""<< folderPng1[0]  << "\" alt=\"The benchmark set SVs count statistics\" style=\"width: 60%;\">\n";
		htmlFile << "</div>\n";*/
		htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
		htmlFile <<"<h4 style=\" margin-left: 40px;\">(a)</h4>\n";
		htmlFile << "  <img src=\""<< folderPng1[0] << "\" alt=\" The benchmark set SVs count statistics\" style=\"width: 50%;\">\n";
		htmlFile <<"<h4 style=\" margin-left: 10px;\">(b)</h4>\n";
		htmlFile << "  <img src=\""<< folderPng1[1]  << "\" alt=\" The benchmark set filtered SVs count statistics\" style=\"width: 50%;\">\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" The quantity distribution of the benchmark set</h4>\n";
		htmlFile <<"</div>\n";
		++figure_num;
		htmlFile <<"The figure shows the distribution of SVs counts of the benchmark set, where (a) represents the result statistics before filtering SVs, and (b) shows the result statistics after filtering SVs.\n";
		htmlFile <<"<br>\n";
		//htmlFile << "The user-called set：<br>";
		for(size_t i=1; i<allmetrics.size(); i++){
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 5px;\">("<< i+1 <<") Statistics of the count of different SV lengths in the user-called set ("<< allmetrics[i][0] <<"): " << "</h4>";
			htmlFile << "The SV reference region size statistics before filtering for user-called set ("<< allmetrics[i][0]<<"):"<< "Total SV count：" << before_filtered[1]<<"<br>";
			htmlFile << "The SV reference region size statistics after filtering for user-called set ("<< allmetrics[i][0]<<"):"<< "Total SV count：" << allmetrics[i][2]<<"<br>";
			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
			htmlFile <<"<h4 style=\" margin-left: 40px;\">(a)</h4>\n";
			htmlFile << "  <img src=\""<< folderPng1[2*(i-1)+2] << "\" alt=\" The user-called set SVs count statistics\" style=\"width: 50%;\">\n";
			htmlFile <<"<h4 style=\" margin-left: 10px;\">(b)</h4>\n";
			htmlFile << "  <img src=\""<< folderPng1[2*(i-1)+3]  << "\" alt=\" The user-called set filtered SVs count statistics\" style=\"width: 50%;\">\n";
			htmlFile << "</div>\n";
			htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
			htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure "<< figure_num  <<" The quantity distribution of the user-called set</h4>\n";
			htmlFile <<"</div>\n";
			++figure_num;
			htmlFile << "The figures show the distribution of SVs counts of the user-called set, where (a) represents the result statistics before filtering SVs, and (b) shows the result statistics after filtering SVs.\n";
		}
		htmlFile << "</div>\n";

		htmlFile << "<br>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<h2>More information</h2>\n";
		htmlFile << "</div>\n";
		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
		htmlFile << "<ul>\n";
		htmlFile << "<li>For more detailed benchmarking results, please refer to the generated result information in the respective folders.</li>\n";
		htmlFile << "<li>For more detailed experiment information, please refer to the github repositories: " << " <a href=\"https://github.com/zhuxiao/asvbm\" target=\"_blank\">asvbm</a>" <<  " and <a href=\"https://github.com/zhuxiao/asvbm-experiments\" target=\"_blank\">asvbm-experiments</a>.</li>\n";
		htmlFile << "<li>For any problems, comments, or suggestions, please contact xzhu@ytu.edu.cn without hesitation. Thank you very much!</li>\n";
		htmlFile << "</ul>\n";
		htmlFile << "-------------------------- This is the end of the Benchmarking Reports. --------------------------\n";
		htmlFile << "</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h2>2. Benchmarking results</h2>";
//		htmlFile <<"</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "Categorizing the benchmarking metrics into two categories after filtering long SV regions: one highlighting performance metrics including recall, precision, and F1 score, and the other presenting benchmark results, encompassing TP_bench, TP_user, FP, FN.   Visualizing these metrics through bar charts offers a more intuitive representation of the assessment results for the detection method.";
//		htmlFile << "</div>";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">The benchmarking results of the user-called set are as follows:</h4>";
//		htmlFile << "</div>\n";
//		// Add h2 style for the first title
//		htmlFile << "<h3>Table 1 Structural variant Detection Method Performance Benchmarking</h3>\n";
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
//		htmlFile << "The table 1 shows the benchmarking results of the structural variant detection method. Where SVs_bench represents the quantity of structural variants (SVs) in the benchmark set, SV_user represents the quantity of SVs in the called set, and SV_user_filtered represents the quantity of SVs after filtering out large SVs. TP stands for True Positive, indicating the count of correctly identified targets or events. FP stands for False Positive, representing the count of falsely identified targets or events. FN represents False Negative, referring to the count of targets or events that were missed or not identified correctly.\n";
//		htmlFile << "</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">The benchmarking results of two categorizes of metrics are shown in the figure:</h4>\n";
//		htmlFile << "</div>\n";
////		htmlFile << "<h3>Figure 1. Benckmarking results of the user-call set</h3>\n";
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
//		htmlFile << "  <img src=\""<< figuresFilePath  << "/benchmarking_result.png\" alt=\"Benchmark results between different tools\">\n";
//		htmlFile << "  <img src=\""<< figuresFilePath  << "/result_classification.png\" alt=\"Performance comparsion between different tools\">\n";
//		htmlFile << "</div>\n";
//		htmlFile << "<div style=\"text-align: center; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h4 style=\"margin-top: 5px; margin-bottom: 0px;\">Figure 1. Benchmarking results of the user-call set</h4>\n";
//		htmlFile <<"</div>\n";

//		// Add h2 style for the second title
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "<h2>3. Statistical results of deviations for overlapping variants</h2>\n";
//		htmlFile <<"</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "For variants that overlap between the user-called set and the benchmark set, the deviations between them are quantified by calculating the center distance and the region size ratio of the overlapping variants.<br>";
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
//		htmlFile << "Calculating the region size ratio for two overlapping variants based on the length of SVs, the closer the ratio is to 1, the smaller the deviation, indicating a more precise and accurate identification result. Statistics results for nine size regions are presented in Table 3:";
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
//		htmlFile << "<h2>4. Benchmarking results for metrics of different SV size regions</h2>\n";
//		htmlFile <<"</div>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "The SV identification results typically contain variants of various sizes, and categorize these variants into different size ranges could be used to explore the identification results more detailed in a fine-grained manner, and could provide new insights into the sensitivity of SV callers to variants of different sizes. Detailed benchmarking results are presented in the table as follows：\n";
//		htmlFile << "</div>\n";
//		//The metric benchmarking results of cutesv in different SV regions
//		for (size_t tableIndex = 0; tableIndex < allregionmetrics.size(); ++tableIndex) {
////			("<<tableIndex+1 <<")
//			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//			htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">("<<tableIndex+1 <<") The user-called set ("<< alltoolnames[tableIndex] <<") of basic metrics results statistics</h3>";
//			htmlFile << "</div>\n";
//			// Add h2 style for the title of each table
//			htmlFile << "<h3> Table "<<(size_t)1 + tableIndex << " The metric benchmarking results of " << alltoolnames[tableIndex] << " in different SV regions</h3>\n";
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
//			htmlFile << "<h4>Benchmarking results for metrics of different SV size regions show as follow with figures:</h4>";
//			htmlFile << "</div>\n";
//			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
//				htmlFile << "  <img src=\""<< folderPng4[2*tableIndex] << "\" alt=\" calculation_result.png\" style=\"width: 30%;\">\n";
//				htmlFile << "  <img src=\""<< folderPng4[2*tableIndex+1]  << "\" alt=\" quantitative_Statistics.png\" style=\"width: 30%;\">\n";
//			htmlFile << "</div>\n";
//		}
//
//		if(folderPng5.size()>0){
//			htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//			htmlFile << "<h3 style=\"margin-top: 5px; margin-bottom: 5px; text-align: left;\">("<< allregionmetrics.size() +1 <<") Benchmarking results for metrics of different SV size regions with different methods</h3>";
//			htmlFile <<"variants are categorized into seven size regions and metrics are computed for comprehensive benchmarking for different detection methods within each region. The benchmarking results are as follows:\n";
//			htmlFile << "</div>\n";
//			htmlFile <<"<div style=\" display: flex; justify-content: center; text-align:center; margin: 0 auto;\">\n";
//				htmlFile << "  <img src=\""<< folderPng5[0] << "\" alt=\"different_range.png\" style=\"width: 50%;\">\n";
//			htmlFile << "</div>\n";
//		}
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile << "For more detailed benchmarking results, please refer to the generated result information in the respective folders.<br>\n";
//		htmlFile << "</div>\n";
//		htmlFile << "<h2>5. Benchmarking Results with Bar Chart </h2>\n";
//		htmlFile << "<div style=\"text-align: left; margin: 0 auto; width: 50%;\">\n";
//		htmlFile <<"Categorizing the evaluation metrics into two categories: one highlighting performance metrics including recall, precision, and F1 score, and the other presenting benchmark results, encompassing TP_bench, TP_user, FP, FN.   Visualizing these metrics through bar charts offers a more intuitive representation of the assessment results for the detection method. The benchmarking results are shown in the figure:\n";
//		htmlFile << "</div>\n";
//		htmlFile << "<h3>Comprehensive benchmarking results</h3>\n";
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
//		htmlFile << "  <img src=\""<< figuresFilePath  << "/benchmarking_result.png\" alt=\"Benchmark results between different tools\">\n";
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
