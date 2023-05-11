#include "iostream"
#include "fstream"
#include "vector"
#include "PYXAID/src_cpp/json.hpp"
using json = nlohmann::json;
using namespace std;
// read params from params.json
int main()
{
    json j;
    ifstream i("params.json");
    i >> j;
    vector<vector<int>> a = j["iconds"].get<vector<vector<int>>>();
    cout << a[2][0] << endl;
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < a[i].size(); j++)
        {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}