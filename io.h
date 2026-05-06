//
// Created by Song Jialong on 24-5-7.
//

#ifndef IO_H
#define IO_H

#include <fstream>
#include <array>
#include<iomanip>

//#include <eigen3/Eigen/Dense>

// 特化的printElelements函数，用于输出std::array
template <typename T, size_t N>
void printElelements(std::ofstream& out, const char* separator, const std::array<T, N>& elements)
{
    out << elements[0];
    for (size_t i = 1; i < elements.size(); ++i)
        out << separator << elements[i];
}

// 特化的printElelements函数，用于输出Eigen::Vector
//template <typename T, int N>
//void printElelements(std::ofstream& out, const char* separator, const Eigen::Matrix<T, N, 1>& elements)
//{
//    out << elements[0];
//    for (size_t i = 1; i < N; ++i)
//        out << separator << elements[i];
//}

// 模板函数printElelements
template <typename T>
void printElelements(std::ofstream& out, const char* separator, const T& elements)
{
    out << elements;
}

// 辅助函数，用于打印容器的某行（递归终止函数）
template <typename Container>
void printLine(std::ofstream& out, const char* separator, size_t line, const Container& container)
{
    // 调用对应的printElelements函数
    printElelements(out, separator, container[line]);
}

// 辅助函数，用于打印容器的某行（递归函数）
template <typename First, typename... Rest>
void printLine(std::ofstream& out, const char* separator, size_t line,
               const First& firstContainer, const Rest&... restContainers)
{
    printElelements(out, separator, firstContainer[line]);
    out << separator;
    printLine(out, separator, line, restContainers...);
}

// 将若干个容器写入文件，要求每个容器具有相同的size
template <typename First, typename... Containers>
void printToFile(std::ofstream& out, const char* separator, const First& first, const Containers&... containers)
{
    size_t lineNum = first.size();
    for (size_t i = 0; i < lineNum; ++i)
    {
        printLine(out, separator, i, first, containers...);
        out << "\n";
    }
}

inline void printCTOC13(std::fstream& fout, int lineID, int launchID, int satID, int event, double time,
                        int targetN1, int targetN2, int targetN3, const double* rv, const double* dv)
{
    fout << std::resetiosflags(std::ios::floatfield | std::ios::adjustfield);
    fout << lineID << " " << launchID << " " << satID << " " << event << " ";
    fout << std::setw(24) << std::setiosflags(std::ios::scientific) << std::setprecision(14) << time << " ";
    fout << std::resetiosflags(std::ios::floatfield | std::ios::adjustfield);
    fout << targetN1 << " " << targetN2 << " " << targetN3 << " ";
    for (int i = 0; i < 6; i++)
    {
        fout << std::setw(24) << std::setiosflags(std::ios::scientific) << std::setprecision(14) << rv[i] / 1e3 << " ";
    }
    fout << std::setw(24) << std::setiosflags(std::ios::scientific) << std::setprecision(14) << dv[0] / 1e3 << " " << dv
        [1] / 1e3 << " " << dv[2] / 1e3 << std::endl;
};


#endif //IO_H
