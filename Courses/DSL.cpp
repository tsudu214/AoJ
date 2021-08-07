#include <iostream>
#include <iomanip>
#include <algorithm>
#include <climits> 
#include <vector>
#include <deque>
#include <list>
#include <stack>
#include <queue>
#include <set>
#include <map>
#include <bitset>
#include <string>
#include <cstring>
#include <cmath>
#include <array>
#include <numeric>

using namespace std;

class DisjointSet
{
public:
    DisjointSet(){}
    DisjointSet(int size){
        rank.resize(size, 0);
        p.resize(size, 0);
        for (int i = 0; i < size; i++) {
            makeSet(i);
        }
    }

    void makeSet(int x)
    {
        p[x] = x;
        rank[x] = 0;
    }

    void unite(int x, int y)
    {
        x = findSet(x);
        y = findSet(y);
        if (rank[x] > rank[y]) {
            p[y] = x;
        }
        else {
            p[x] = y;
            if (rank[x] == rank[y]) {
                rank[y]++;
            }
        }
    }

    int findSet(int x) 
    {
        if (p[x] != x) {
            p[x] = findSet(p[x]);
        }
        return p[x];
    }

private:
    vector<int> rank;
    vector<int> p;

};


int main()
{
    int n, q;
    cin >> n >> q;

    DisjointSet ds(n);

    for (int i = 0; i < q; i++) {
        int com, x, y;
        cin >> com >> x >> y;
        if (com == 0) {
            ds.unite(x, y);
        }
        else if (com == 1) {
            if (ds.findSet(x) == ds.findSet(y)) {
                cout << 1 << endl;
            } else {
                cout << 0 << endl;
            }
        }
    }

    return 0;
}

