#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <stdexcept>

using namespace std;
int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: ./a.out <in> <out>" << endl;
        return 1;
    }
    if (std::ifstream(argv[2])) {
      std::cerr << "target file already exists" << std::endl;
      return 1;
    }
    ifstream ifs(argv[1]);
    ofstream ofs(argv[2], std::ios::binary);

    std::cout << "Reading matrix \"" << argv[1] << "\"" << std::endl;
    // Skip %% MatrixMarket header and any comments.
    char line[256];
    char percent[1];
    percent[0] = '%';
    ifs.getline(line, 256);
    while (strncmp(line,percent,1) == 0) {
      ifs.getline(line, 256);
    }

    int m, n, nnz;
    int numConverted = sscanf(line,"%d %d %d",&m,&n,&nnz);

    if (numConverted != 3) {
      std::ostringstream errStr;
      errStr << "Error reading matrix dimensions.  Expected 3 integers, found " << numConverted;
      throw(std::runtime_error(errStr.str()));
    }
    ofs.write(reinterpret_cast<char*>(&m),   sizeof(m));
    ofs.write(reinterpret_cast<char*>(&n),   sizeof(n));
    ofs.write(reinterpret_cast<char*>(&nnz), sizeof(nnz));

    int    row = -99, i, j;
    double v;
    vector<int>    inds;
    vector<double> vals;

    ifs >> i >> j >> v;
    i--; j--;
    while (row != i && ifs.good()) {
        row = i;

        inds.resize(0);
        vals.resize(0);

        do {
            inds.push_back(j);
            vals.push_back(v);

            ifs >> i >> j >> v;
            i--; j--;
        } while (row == i && ifs.good());


        int rownnz = inds.size();
        ofs.write(reinterpret_cast<char*>(&row),    sizeof(row));
        ofs.write(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
        for (int k = 0; k < rownnz; k++) ofs.write(reinterpret_cast<char*>(&inds[0] + k), sizeof(inds[k]));
        for (int k = 0; k < rownnz; k++) ofs.write(reinterpret_cast<char*>(&vals[0] + k), sizeof(vals[k]));

    }

    return 0;
}
