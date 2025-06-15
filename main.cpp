
#include "multiFSPBWT.h"




void printHelp(const char* programName) {
    std::cout << "用法: " << programName << " [-i input] [-o output] [-B block_size] [-f F_value] [-l length]\n"
              << "  -i <file>  输入文件 (默认: out)\n"
              << "  -o <file>  输出文件 (默认: <input>.out)\n"
              << "  -B <int>   音节大小 B (默认: 64)\n"
              << "  -f <int>   F 值 (默认: 1)\n"
              << "  -l <int>   查询长度 (默认: 1800)\n"
              << "示例: " << programName << " -i data.txt -B 128 -l 2000\n";
}

// 验证文件路径
bool validateFiles(const std::string& inputFile, const std::string& outputFile) {
    std::ifstream in(inputFile);
    if (!in.good()) {
        std::cerr << "错误: 输入文件 '" << inputFile << "' 不存在\n";
        return false;
    }
    in.close();
    std::ofstream out(outputFile, std::ios::app);
    if (!out.good()) {
        std::cerr << "错误: 无法写入输出文件 '" << outputFile << "'\n";
        return false;
    }
    out.close();
    return true;
}

int main(int argc, char* argv[]) {
    std::string inputFile = "out";    // 默认输入文件名
    std::string outputFile;           // 输出文件动态生成
    int B = 64;                       // 默认 B 值
    int F = 1;                        // 默认 F 值
    int queryLength = 1800;           // 默认查询长度

    // 解析命令行参数
    int opt;
    while ((opt = getopt(argc, argv, "i:I:o:O:b:B:f:F:l:L:h:H")) != -1)
    {
        try {
            switch (opt) {
            case 'i':
            case 'I':
                inputFile = optarg;
                break;
            case 'o':
            case 'O':
                outputFile = optarg;
                break;
            case 'b':
            case 'B':
                B = std::stoi(optarg);
                break;
            case 'f':
            case 'F':
                F = std::stoi(optarg);
                break;
            case 'l':
            case 'L':
                queryLength = std::stoi(optarg);
                break;
            case 'h':
            case 'H':
                printHelp(argv[0]);
                return 0;
            default:
                std::cerr << "错误: 未知选项，使用 -h 查看用法\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "错误: 无效的参数值 '" << optarg << "'\n";
            return 1;
        }
    }
    // 自动生成输出文件名（若未指定）
    if (outputFile.empty()) {
        outputFile = inputFile + ".out";
    }

    // 验证参数
    if (B <= 0 || F < 0 || queryLength <= 0) {
        std::cerr << "错误: B、F 和查询长度必须为正数\n";
        return 1;
    }
    if (!validateFiles(inputFile, outputFile)) {
        return 1;
    }

    // 输出参数信息
    std::cout << "参数:\n"
              << "输入: " << inputFile << "\n"
              << "输出: " << outputFile << "\n"
              << "B: " << B << "\n"
              << "F: " << F << "\n"
              << "查询长度: " << queryLength << "\n";

    // 执行程序逻辑
    multiFSPBWT<unsigned long long> CRY;
    CRY.F = F;
    CRY.B = B;
    CRY.T = pow(2, F);
    CRY.minSiteL = B * 2 - 1;

    int a = CRY.readMacsPanel(inputFile);
    std::cout << "读取面板: " << a << "\n";
    if (a != 0) return a;

    int b = CRY.makeFuzzyPanelGlobally();
    std::cout << "生成模糊面板: " << b << "\n";
    if (b != 0) return b;

    int c = CRY.inPanelLongMatchQuery(queryLength, outputFile);
    std::cout << "查询完成: " << c << "\n";
    return c;
}