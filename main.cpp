
#include "multiFSPBWT.h"

void printHelp()
{
	std::cout << "Usage: program [options]\n";
	std::cout << "Options:\n";
	std::cout << "  -h, -H             Print this help message and exit\n";
	std::cout << "  -B, -b <value>     Set the value of B (default: 64)\n";
	std::cout << "  -F, -f <value>     Set the value of F (default: 2)\n";
	std::cout
			<< "  -i, -I <file>      Specify the panel file (default: panel.vcf)\n";
	std::cout
			<< "  -o, -O <file>      Specify the output file (default: panel_B_F_L.txt)\n";
	std::cout << "  -L, -l <value>     Set the value of L (default: 500)\n";
	std::cout
			<< "  -m, -M <mode>      in for in-panel query; out for out-panel query\n";
	std::cout << "  -q, -Q <file>      Specify the query file(default: panel.vcf)\n";
	std::cout << "  -e, -even          Set the even(position-based) mode (default: count-based)\n";
	std::cout << std::endl;
}

/*
int main(int argc, char *argv[])
{
	int B = 64, F = 2, L = 500;
	string panelFile = "panel.vcf";
	string outputFile = "";
	string mode = "in";
	string queryFile="query.vcf";
	bool save = false;
	string saveFile = "";
	bool load = false;
	string loadFile = "";
	bool even = false;

	if (argc == 1)
	{
		printHelp();
		return 0;
	}

	for (int i = 1; i < argc; i++)
	{
		string arg = argv[i];
		if (arg == "-h" || arg == "-H")
		{
			printHelp(); // 打印帮助信息
			exit(0); // 退出程序
		}
		if (arg == "-B" || arg == "-b")
		{
			if (i + 1 < argc)
			{
				B = atoi(argv[i + 1]);
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-F" || arg == "-f")
		{
			if (i + 1 < argc)
			{
				F = atoi(argv[i + 1]);
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-i" || arg == "-I")
		{
			if (i + 1 < argc)
			{
				panelFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-o" || arg == "-O")
		{
			if (i + 1 < argc)
			{
				outputFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-L" || arg == "-l")
		{
			if (i + 1 < argc)
			{
				L = atoi(argv[i + 1]);
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-m" || arg == "-M")
		{
			if (i + 1 < argc)
			{
				mode = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-q" || arg == "-Q")
		{
			if (i + 1 < argc)
			{
				queryFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}

		else if (arg == "-save" || arg == "-SAVE")
		{
			save = true;
			if (i + 1 < argc)
			{
				saveFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-load" || arg == "-LOAD")
		{
			load = true;
			if (i + 1 < argc)
			{
				loadFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-e" || arg == "-even")
		{
			even = true; // 如果出现 -e 或 -even 选项，则将 even 设置为 true
		}
	}



	if (F<=0 || F>4) {
		cout << "wrong fuzzy size! must be 1/2/3/4" << endl;
	}
	if (outputFile=="")
	{
		if (mode == "out") {
			outputFile = "FSPBWT_outPanel_" + panelFile + "_" + queryFile+ "_" + std::to_string(B)
					+ "_" + std::to_string(F) + "_" + std::to_string(L)
					+ ".txt";
		}
		else if (mode == "in") {
			outputFile = "FSPBWT_inPanel_" + panelFile + "_" + std::to_string(B)
					+ "_" + std::to_string(F) + "_" + std::to_string(L)
					+ ".txt";
		}

	}
	string informationFile="Inf_FSPBWT_";
	if (mode == "out") {
		informationFile += panelFile + "_outPanelQuery_" + std::to_string(B) + "_"
				+ std::to_string(F) + "_" + std::to_string(L) + ".txt";
	}
	else if (mode == "in") {
		informationFile += panelFile + "_inPanelQuery_" + std::to_string(B) + "_"
					+ std::to_string(F) + "_" + std::to_string(L) + ".txt";

	}
	else {
		std::cout << "wrong mode! must be in / out" << endl;
	}

	// 输出参数
	cout << "Parameters:" << endl;
	cout << "  B: " << B << endl;
	cout << "  F: " << F << endl;
	cout << "  L: " << L << endl;
	cout << "  Input file: " << panelFile << endl;
	cout << "  Output file: " << outputFile << endl;
	cout << "  Mode: " << mode << endl;
	cout << "  Query file: " << queryFile << endl;
	cout << "  Save: " << (save ? "true" : "false") << endl;
	cout << "  Save file: " << saveFile << endl;
	cout << "  Load: " << (load ? "true" : "false") << endl;
	cout << "  Load file: " << loadFile << endl;
	cout << "  Even: " << (even ? "true" : "false") << endl;

	if (B==64) {
		multiFSPBWT<unsigned long long> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值
		int a = CRY.readVCF(panelFile);
		// int a = CRY.readTXT("sites.txt");
			std::cout << "read panel file done: " << a << endl;
			//1846144
		int b;
		if (even==true)
		{
			b = CRY.makeFuzzyPanelEvenly();
		}
		else
		{
			b = CRY.makeFuzzyPanelGlobally();

		}
			//2139264
			std::cout << "make fuzzy panel  done: " << b << endl;
			if (mode == "out") {
				int c = CRY.readQueryVCF(queryFile);
				cout << "read query done: " << c << endl;

				// 2146176
				int d = CRY.outPanelLongMatchQuery(L, outputFile, even);
				cout << "out-panel query done: " << d << endl;
				// 2146688
				CRY.outputInformationToFile(informationFile, "out");
			}
			else if (mode == "in") {
				int c = CRY.inPanelLongMatchQuery(L, outputFile);
				std::cout << "in-panel query done: " << c << endl;
				CRY.outputInformationToFile(informationFile, "in");
			}
			return 0;
	}
	else if (B==128) {
		multiFSPBWT<unsigned __uint128_t> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值
		int a = CRY.readVCF(panelFile);

			std::cout << "read panel file done: " << a << endl;
		int b;
		if (even==true)
		{
			b = CRY.makeFuzzyPanelEvenly();
		}
		else
		{
			b = CRY.makeFuzzyPanelGlobally();

		}
			std::cout << "make fuzzy panel  done: " << b << endl;

			if (mode == "out") {
				int c = CRY.readQueryVCF(queryFile);
				cout << "read query done: " << c << endl;
				int d = CRY.outPanelLongMatchQuery(L, outputFile, even);
				cout << "out-panel query done: " << d << endl;
				CRY.outputInformationToFile(informationFile, "out");
			}
			else if (mode == "in") {
				int c = CRY.inPanelLongMatchQuery(L, outputFile);
				std::cout << "in-panel query done: " << c << endl;
				CRY.outputInformationToFile(informationFile, "in");
			}
			return 0;

	}
	else {
		std::cout << "wrong Syllable size! s must be 64 / 128" << endl;
	}
	return 0;
}


*/

void printMultiFSPBWTStats(const multiFSPBWT<unsigned long long>& cry) {
    // 输出基本配置参数
    std::cout << "===== 配置参数 =====" << std::endl;
    std::cout << "B: " << cry.B << std::endl;
    std::cout << "F: " << cry.F << std::endl;
    std::cout << "T: " << cry.T << std::endl;
    std::cout << "M: " << cry.M << std::endl;
    std::cout << "N: " << cry.N << std::endl;
    std::cout << "n: " << cry.n << std::endl;
    std::cout << "minSiteL: " << cry.minSiteL << std::endl;

    // 输出时间统计
    std::cout << "\n===== 时间统计(秒) =====" << std::endl;
    std::cout << "readPanelTime: " << cry.readPanelTime << std::endl;
    std::cout << "makeFuzzyPanelTime: " << cry.makeFuzzyPanelTime << std::endl;
    std::cout << "inPanelQueryTime: " << cry.inPanelQueryTime << std::endl;
    std::cout << "inPanelFilteringTime: " << cry.inPanelFilteringTime << std::endl;
    std::cout << "inPanelIdentificationTime: " << cry.inPanelIdentificationTime << std::endl;
    std::cout << "readQueryTime: " << cry.readQueryTime << std::endl;
    std::cout << "outPanelQueryTime: " << cry.outPanelQueryTime << std::endl;
    std::cout << "outPanelFilteringTime: " << cry.outPanelFilteringTime << std::endl;
    std::cout << "outPanelIdentificationTime: " << cry.outPanelIdentificationTime << std::endl;

    // 输出计数统计
    std::cout << "\n===== 计数统计 =====" << std::endl;
    std::cout << "inPanelAlternativeNum: " << cry.inPanelAlternativeNum << std::endl;
    std::cout << "inPanelMatchNum: " << cry.inPanelMatchNum << std::endl;
    std::cout << "outPanelAlternativeNum: " << cry.outPanelAlternativeNum << std::endl;
    std::cout << "outPanelMatchNum: " << cry.outPanelMatchNum << std::endl;
    std::cout << "alternativeSyllableNum: " << cry.alternativeSyllableNum << std::endl;
    std::cout << "matchLen: " << cry.matchLen << std::endl;
    std::cout << "panelMultiSyllableNum: " << cry.panelMultiSyllableNum << std::endl;

    // 计算并输出比率
    std::cout << "\n===== 计算比率 =====" << std::endl;
    long long product_M_n = static_cast<long long>(cry.M) * cry.n;
    std::cout << "M * n: " << product_M_n << std::endl;

    if (product_M_n > 0) {
        double ratio = static_cast<double>(cry.panelMultiSyllableNum) / product_M_n;
        std::cout << "panelMultiSyllableNum / (M * n): " << ratio << std::endl;
    } else {
        std::cout << "panelMultiSyllableNum / (M * n): 无法计算 (M*n为0)" << std::endl;
    }
	// 输出 panelCount 数组
	std::cout << "\n===== Panel Count 统计 =====" << std::endl;
	for (int i = 0; i < 10; ++i) {
		std::cout << i << "的数量: " << cry.panelCount[i] << std::endl;
	}
}

template<typename Syllable>
void printMemoryUsage(const multiFSPBWT<Syllable>& instance) {
    std::cout << "\n===== 内存占用统计 =====" << std::endl;

    // 计算 X 的内存占用
    size_t xSize = 0;
    for (const auto& row : instance.X) {
        xSize += row.capacity() * sizeof(Syllable);
    }
    std::cout << "X: " << xSize / (1024.0 * 1024.0) << " MB (" << xSize << " bytes)" << std::endl;

    // 计算 fuzzyX 的内存占用
    size_t fuzzyXSize = 0;
    for (const auto& row : instance.fuzzyX) {
        fuzzyXSize += row.capacity() * sizeof(uint32_t);
    }
    std::cout << "fuzzyX: " << fuzzyXSize / (1024.0 * 1024.0) << " MB (" << fuzzyXSize << " bytes)" << std::endl;

    // 计算 array 的内存占用
    size_t arraySize = 0;
    for (const auto& row : instance.array) {
        arraySize += row.capacity() * sizeof(int);
    }
    std::cout << "array: " << arraySize / (1024.0 * 1024.0) << " MB (" << arraySize << " bytes)" << std::endl;

    // 计算 divergence 的内存占用
    size_t divergenceSize = 0;
    for (const auto& row : instance.divergence) {
        divergenceSize += row.capacity() * sizeof(int);
    }
    std::cout << "divergence: " << divergenceSize / (1024.0 * 1024.0) << " MB (" << divergenceSize << " bytes)" << std::endl;

    // 计算 panelMultiSyllable 的内存占用
    size_t panelMultiSyllableSize = 0;
    for (const auto& row : instance.panelMultiSyllable) {
        panelMultiSyllableSize += row.capacity() * sizeof(bool);
    }
    std::cout << "panelMultiSyllable: " << panelMultiSyllableSize / (1024.0 * 1024.0) << " MB (" << panelMultiSyllableSize << " bytes)" << std::endl;

    // 计算 queryMultiSyllable 的内存占用
    size_t queryMultiSyllableSize = 0;
    for (const auto& row : instance.queryMultiSyllable) {
        queryMultiSyllableSize += row.capacity() * sizeof(bool);
    }
    std::cout << "queryMultiSyllable: " << queryMultiSyllableSize / (1024.0 * 1024.0) << " MB (" << queryMultiSyllableSize << " bytes)" << std::endl;

    // 计算 u 的内存占用
    size_t uSize = instance.n * instance.M * instance.T * sizeof(int);
    std::cout << "u: " << uSize / (1024.0 * 1024.0) << " MB (" << uSize << " bytes)" << std::endl;

    // 计算 panelMultiMaps 的内存占用 (std::vector<std::vector<Uint4Array>>)
    size_t panelMultiMapsSize = 0;
    for (const auto& row : instance.panelMultiMaps) {
        panelMultiMapsSize += sizeof(std::vector<Uint4Array>) + row.capacity() * sizeof(Uint4Array);
        for (const auto& arr : row) {
            panelMultiMapsSize += arr.get_data().capacity() * sizeof(uint8_t);
        }
    }
    std::cout << "panelMultiMaps: " << panelMultiMapsSize / (1024.0 * 1024.0) << " MB (" << panelMultiMapsSize << " bytes)" << std::endl;

    // 计算 queryMultiMaps 的内存占用 (std::vector<std::unordered_map<int, std::string>>)
    size_t queryMultiMapsSize = 0;
    for (const auto& map : instance.queryMultiMaps) {
        queryMultiMapsSize += map.size() * (sizeof(std::pair<int, std::string>) +
                                           map.bucket_count() * sizeof(void*));
        for (const auto& pair : map) {
            queryMultiMapsSize += pair.second.capacity() * sizeof(char);
        }
    }
    std::cout << "queryMultiMaps: " << queryMultiMapsSize / (1024.0 * 1024.0) << " MB (" << queryMultiMapsSize << " bytes)" << std::endl;

    // 计算总内存占用
    size_t totalSize = xSize + fuzzyXSize + arraySize + divergenceSize +
                       panelMultiSyllableSize + queryMultiSyllableSize +
                       uSize + panelMultiMapsSize + queryMultiMapsSize;
    std::cout << "总计: " << totalSize / (1024.0 * 1024.0) << " MB (" << totalSize << " bytes)" << std::endl;
}
int main(int argc, char *argv[])
{
	int B=64,F=1;
	multiFSPBWT<unsigned long long> CRY;
	CRY.F=F;
	CRY.B = B;
	CRY.T = pow(2, F); // 计算T的值，即2的F次方
	CRY.minSiteL = B * 2 - 1;
	// 构造函数用于初始化B、F、T和minSiteL的值
	string panel = "out";
	int a = CRY.readMacsPanel(panel);
	std::cout << "read panel file done: " << a << endl;
	int b = CRY.makeFuzzyPanelGlobally();
	std::cout << "make panel file done: " << b << endl;
	int c = CRY.inPanelLongMatchQuery(1800,"outFile");
	std::cout << "query done: " << c << endl;
	//
	printMultiFSPBWTStats(CRY);
	printMemoryUsage(CRY);
	return 0;
}
