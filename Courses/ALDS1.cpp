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

#ifdef ALDS1_12_C

const int INF = 1 << 21;

struct edge {
    int to;
    int cost;
};

vector<vector<edge>> g;

vector<int> dijkstra(int s, int n)
{
    vector<int> p(n, -1);
    vector<int> d(n, INF);

    auto c = [&](int l, int r) { return d[l] > d[r]; };
    priority_queue<int, vector<int>, decltype(c) > PQ(c);

    d[s] = 0;
    PQ.push(s);

    while(!PQ.empty()) {
        int u = PQ.top(); PQ.pop();

        for ( auto& e : g[u]) {
            if ( d[u] + e.cost < d[e.to] ) {
                d[e.to] = d[u] + e.cost;
                p[e.to] = u;
                PQ.push(e.to);
            }
        }
    }

    return d;
}

int main()
{
    int n;
    cin >> n;
    g.resize(n);

    for (int i = 0; i < n; i++) {
        int u, k;
        cin >> u >> k;
        for (int j = 0; j < k; j++) {
            int v, c;
            cin >> v >> c;
            edge e = { v, c };
            g[u].push_back(e);
        }
    }

    vector<int> d = dijkstra(0, n);

    for ( int i = 0; i < n; i++ ) {
        cout << i << " " << d[i] << endl;
    }

    return 0;
}

#endif

#ifdef ALDS1_12_B

int g[200][200];
const int INF = 1 << 21;
enum { WHITE=0, GRAY=1, BLACK=2 };

void dijkstra(int s, int n, vector<int>& d)
{
    vector<int> color(n, WHITE);
    vector<int> p(n, -1);

    d[s] = 0;

    while(true) {
        int u = -1;
        int min = INF;
        for (int i = 0; i < n; i++) {
            if ( color[i] != BLACK && d[i] < min ) {
                u = i;
                min = d[i];
            }
        }
        if (u < 0) {
            break;
        }

        color[u] = BLACK;
        for (int v = 0; v < n; v++) {
            if ( color[v] != BLACK && g[u][v] < INF ) {
                if ( d[u] + g[u][v] < d[v] ) {
                    d[v] = d[u] + g[u][v];
                    p[v] = u;
                    color[v] = GRAY;
                }
            }
        }
    }
}

int main()
{
    int n;
    cin >> n;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            g[i][j] = INF;
        }
    }

    for (int i = 0; i < n; i++) {
        int u, k;
        cin >> u >> k;
        for (int j = 0; j < k; j++) {
            int v, c;
            cin >> v >> c;
            g[u][v] = c;
        }
    }

    vector<int> d(n, INF);

    dijkstra(0, n, d);

    for ( int i = 0; i < n; i++ ) {
        cout << i << " " << d[i] << endl;
    }

    return 0;
}

#endif

#ifdef ALDS1_12_A

int g[200][200];
const int INF = 1 << 21;
enum { WHITE=0, GRAY=1, BLACK=2 };

int prim(int s, int n)
{
    vector<int> color(n, WHITE);
    vector<int> p(n, -1);
    vector<int> d(n, INF);

    d[s] = 0;

    while(true) {
        int v = -1;
        int min = INF;
        for (int i = 0; i < n; i++) {
            if ( color[i] != BLACK && d[i] < min ) {
                v = i;
                min = d[i];
            }
        }
        if (v < 0) {
            break;
        }

        color[v] = BLACK;
        for (int w = 0; w < n; w++) {
            if ( color[w] != BLACK && g[v][w] < INF ) {
                if ( g[v][w] < d[w] ) {
                    d[w] = g[v][w];
                    p[w] = v;
                    color[w] = GRAY;
                }
            }
        }
    }

    int sum = 0;
    for (int i = 0; i < n; i++) {
        if ( p[i] >= 0 ) {
            sum += g[i][p[i]];
        }
    }
    return sum;
}

int main()
{
    int n;
    cin >> n;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            g[i][j] = INF;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int w;
            cin >> w;
            if (w >= 0) g[i][j] = w;
        }
    }

    int sum = prim(0, n);

    cout << sum << endl;

    return 0;
}

#endif

#ifdef ALDS1_11_D

void dfs(const vector<vector<int>>& g, int s, vector<int>& color)
{
    int dist = 0;
    stack<int> S;
    S.push(s);
    color[s] = s;
    while (!S.empty()) {
        int v = S.top();
        int n = -1;
        for ( auto& w : g[v] ) {
            if (color[w] == -1) {
                n = w;
                break;
            }
        }
        if (n >= 0) {
            S.push(n);
            color[n] = s;
        }
        else {
            S.pop();
        }
    }
}

int main()
{
    int n, e;
    cin >> n >> e;
    vector<vector<int>> g(n);

    for (int i = 0; i < e; i++) {
        int v, w;
        cin >> v >> w;
        g[v].push_back(w);
        g[w].push_back(v);
    }

    vector<int> color(n, -1);

    for (int i = 0; i < n; i++) {
        if (color[i] == -1) {
            dfs(g, i, color);
        }
    }

    int q;
    cin >> q;
    for (int i = 0; i < q; i++ ) {
        int v, w;
        cin >> v >> w;
        cout << ((color[v] == color[w])? "yes" : "no") << endl;
    }

    return 0;
}

#endif

#ifdef ALDS1_11_C

enum { WHITE=0, BLACK=1 };

void bfs(const vector<vector<int>>& g, int s, vector<int>& color, vector<int>& d)
{
    int dist = 0;
    queue<int> Q;
    Q.push(s);
    color[s] = BLACK;
    d[s] = 0;
    while (!Q.empty()) {
        int v = Q.front(); Q.pop();
        for ( auto& w : g[v] ) {
            if (color[w] == WHITE) {
                Q.push(w);
                color[w] = BLACK;
                d[w] = d[v]+1;
            }
        }
    }
}

int main()
{
    int n;
    cin >> n;
    vector<vector<int>> g(n);

    for (int i = 0; i < n; i++) {
        int v, k;
        cin >> v; v--;
        cin >> k;
        for (int j = 0; j < k; j++) {
            int w;
            cin >> w; w--;
            g[v].push_back(w);
        }
        sort(g[v].begin(), g[v].end());
    }

    vector<int> color(n, WHITE);
    vector<int> d(n, -1);

    bfs(g, 0, color, d);

    for (int i = 0; i < n; i++ ) {
        cout << i+1 << " " << d[i] << endl;
    }

    return 0;
}

#endif

#ifdef ALDS1_11_B

enum { WHITE=0, GRAY=1, BLACK=2 };

void dfs(const vector<vector<int>>& g, int s, vector<int>& color, vector<int>& d, vector<int>&f, int& time)
{
    stack<int> S;
    S.push(s);
    color[s] = GRAY;
    d[s] = ++time;
    while (!S.empty()) {
        int v = S.top();
        int n = -1;
        for ( auto& w : g[v] ) {
            if (color[w] == WHITE) {
                n = w;
                break;
            }
        }
        if (n >= 0) {
            S.push(n);
            color[n] = GRAY;
            d[n] = ++time;
        }
        else {
            color[v] = BLACK;
            f[v] = ++time;
            S.pop();
        }
    }
}

int main()
{
    int n;
    cin >> n;
    vector<vector<int>> g(n);

    for (int i = 0; i < n; i++) {
        int v, k;
        cin >> v; v--;
        cin >> k;
        for (int j = 0; j < k; j++) {
            int w;
            cin >> w; w--;
            g[v].push_back(w);
        }
        sort(g[v].begin(), g[v].end());
    }

    vector<int> color(n, WHITE);
    vector<int> d(n, -1);
    vector<int> f(n, -1);
    int time = 0;

    for (int i = 0; i < n; i++) {
        if (color[i] == WHITE) {
            dfs(g, i, color, d, f, time);
        }
    }

    for (int i = 0; i < n; i++ ) {
        cout << i+1 << " " << d[i] << " " << f[i] << endl;
    }

    return 0;
}

#endif

#ifdef ALDS1_11_A

int main()
{
    int n;
    cin >> n;
    vector<vector<bool>> g(n);
    for (int i = 0; i < n; i++) {
        g[i].resize(n, false);
    }

    for (int i = 0; i < n; i++) {
        int v, k;
        cin >> v; v--;
        cin >> k;
        for (int j = 0; j < k; j++) {
            int w;
            cin >> w; w--;
            g[v][w] = true;
        }
    }
    
    for (int i = 0; i < n; i++ ) {
        for (int j = 0; j < n; j++) {
            if (j > 0) cout << " ";
            cout << (g[i][j]? 1 : 0);
        }
        cout << endl;
    }

    return 0;
}

#endif

#ifdef ALDS1_10_D

double cost[1000][1000] = {0};

int main()
{
    
    int n;
    cin >> n;
    vector<double> p(n+1), q(n+1);
    p[0] = 0;
    for (int i = 1; i <= n; i++) {
        cin >> p[i];
    }
    for (int i = 0; i <= n; i++) {
        cin >> q[i];
    }
    
    //int n = 5;
    //vector<double> p = {0, 0.1500, 0.1000, 0.0500, 0.1000, 0.2000};
    //vector<double> q = {0.0500, 0.1000, 0.0500, 0.0500, 0.0500, 0.1000};

    //int n = 7;
    //vector<double> p = { 0, 0.0400, 0.0600, 0.0800, 0.0200, 0.1000, 0.1200, 0.1400};
    //vector<double> q = {0.0600, 0.0600, 0.0600, 0.0600, 0.0500, 0.0500, 0.0500, 0.0500};

    for (int j = 0; j <= n; j++ ) {
        cost[j+1][j] = q[j];
    }
    for (int j = 1; j <= n; j++ ) {
        for (int i = j; i >= 1; i--) {
            cost[i][j] = HUGE_VAL;
            double pqsum = accumulate(p.begin()+i, p.begin()+j+1, (double)0) + accumulate(q.begin()+i-1, q.begin()+j+1, (double)0);
            for (int m = j; m >= i; m--) {
                cost[i][j] = min(cost[i][j], cost[i][m-1] + cost[m+1][j] + pqsum );
            }
        }
    }

    cout << fixed << setprecision(8);
    cout << cost[1][n] << endl;

    return 0;
}

#endif

#ifdef ALDS1_10_C

int dp[1001][1001];

int getlength(const char* x, const char* y)
{
    int m = (int)strlen(x);
    int n = (int)strlen(y);

    bool found = false;
    for (int j = 0; j < n; j++ ) {
        if (x[0] == y[j]) found = true;
        dp[0][j] = found? 1 : 0;
    }
    found = false;
    for (int i = 0; i < m; i++ ) {
        if (y[0] == x[i]) found = true;
        dp[i][0] = found? 1 : 0;
    }

    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++ ) {
            if (x[i] == y[j]) {
                dp[i][j] = dp[i-1][j-1] + 1;
            }
            else {
                dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
            }
        }
    }

    return dp[m-1][n-1];
}

int main()
{
    
    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        string x, y;
        cin >> x;
        cin >> y;
        int L = getlength(x.c_str(), y.c_str());
        cout << L << endl;
    }
    
    return 0;
}

#endif

#ifdef ALDS1_10_B

int dp[101][101];

int main()
{
    int n;
    
    cin >> n;

    vector<int> vn(n+1);
    for (int i = 0; i < n; i++) {
        int r, c;
        cin >> r >> c;
        vn[i] = r;
        if (i == n-1) {
            vn[n] = c;
        }
    }
    
    //n = 6;
    //vector<int> vn = {30, 35, 15, 5, 10, 20, 25};
    //n = 5;
    //vector<int> vn = {1, 34, 44, 13, 30, 1};

    // M[0]M[1]...M[n-1]
    // M[i] : vn[i] x vn[i+1]
    // dp[i][j] = min cost of M[i]..M[j] 0<=i<=n-2, i+1<=j<=n-1

    for (int j = 1; j < n; j++) {
        dp[j][j] = 0;
        for (int i = j-1; i >= 0; i--) {
            int min_cost = 1 << 21;
            for (int m = i; m <= j-1; m++) {
                int cost = dp[i][m] + dp[m+1][j] + vn[i] * vn[m+1] * vn[j+1];
                if (cost < min_cost) min_cost = cost;
            }
            dp[i][j] = min_cost;
        }
    }

    cout << dp[0][n-1] << endl;

    return 0;
}

#endif

#ifdef ALDS1_10_A

int main()
{
    int F[50];
    F[0] = F[1] = 1;

    int n;
    cin >> n;

    for (int i = 2; i <= n; i++ ) {
        F[i] = F[i-1] + F[i-2];
    }

    cout << F[n] << endl;

    return 0;
}

#endif

#ifdef ALDS1_9_C

void maxHeapify(int n, int A[], int i)
{
    int l = 2*i, r = 2*i+1;
    int largest = i;
    if (l <= n && A[l] > A[largest] ) largest = l;
    if (r <= n && A[r] > A[largest] ) largest = r;

    if (largest != i) {
        swap(A[i], A[largest]);
        maxHeapify(n, A, largest);
    }
}

void insert(int* pn, int A[], int k)
{
    int i = ++(*pn);
    A[i] = k;
    while (i > 1 && A[i/2] < A[i]) {
        swap(A[i/2], A[i]);
        i /= 2;
    }
}

int extract(int* pn, int A[])
{
    int n = *pn;
    int M = A[1];
    swap(A[1], A[n]);
    n--;
    maxHeapify(n, A, 1);
    *pn = n;
    return M;
}

int main()
{
    int n = 0;
    vector<int> A(1, 0);
    for (int i = 0; i <= 2000000000; i++) {
        char com[20];
        int k;
        cin >> com;
        if (com[0] == 'i') {
            cin >> k;
            A.push_back(k);
            insert(&n, &A[0], k);
        }
        else if (com[0] == 'e' && com[1] == 'x' ) {
            cout << extract(&n, &A[0]) << endl;
            A.pop_back(); 
        }
        else if (com[0] == 'e' && com[1] == 'n' ) {
            break;
        }
    }

    return 0;
}

#endif

#ifdef ALDS1_9_B

void maxHeapify(int n, int A[], int i)
{
    int l = 2*i, r = 2*i+1;
    int largest = i;
    if (l <= n && A[l] > A[largest] ) largest = l;
    if (r <= n && A[r] > A[largest] ) largest = r;

    if (largest != i) {
        swap(A[i], A[largest]);
        maxHeapify(n, A, largest);
    }
}

void buildMaxHeap(int n, int A[])
{
    for (int i = n/2; i >= 1; i--) {
        maxHeapify(n, A, i);
    }
}

int main()
{
    int n;
    cin >> n;
    vector<int> A(n+1, 0);
    for (int i = 1; i <= n; i++) cin >> A[i];

    buildMaxHeap(n, &A[0]);

    for (int i = 1; i <= n; i++) cout << " " << A[i];
    cout << endl;

    return 0;
}

#endif

#ifdef ALDS1_9_A

int main()
{
    int n;
    cin >> n;
    vector<int> A(n+1, 0);
    for (int i = 1; i <= n; i++) cin >> A[i];

    for (int i = 1; i <= n; i++) {
        int p = i/2, left = 2*i, right = 2*i+1;
        cout << "node " << i << ": key = " << A[i] << ", ";
        if (1 <= p) {
            cout << "parent key = " << A[p] << ", ";
        }
        if (left <= n) {
            cout << "left key = " << A[left] << ", ";
        }
        if (right <= n) {
            cout << "right key = " << A[right] << ", ";
        }
        cout << endl;
    }

    return 0;
}

#endif

#ifdef ALDS1_8_C

struct Node;

Node* NIL;

struct Node {
    Node(int k) : key(k), parent(NIL), left(NIL), right(NIL) {};
    int key;
    Node *parent;
    Node *left, *right;
};

Node* root = NIL;

void insert(int k)
{
    Node* x = root;
    Node* y = NIL;   // parent of x
    Node* z = new Node(k);
    
    while (x != NIL) {
        y = x;
        if (z->key < x->key) {
            x = x->left;
        } else {
            x = x->right;
        }
    }
    z->parent = y;

    if (y == NIL) { // first node
        root = z;
    }
    else {
        if (z->key < y->key) {
            y->left = z;
        } else {
            y->right = z;
        }
    }
}

Node* find(int k)
{
    Node* x = root;
    while (x != NIL) {
        if (k < x->key) {
            x = x->left;
        } 
        else if (x->key < k) {
            x = x->right;
        }
        else {
            return x;
        }
    }
    return NIL;
}

void del(Node* z)
{
    Node* p = z->parent;
    if (p == NIL) {
        root = NIL;
        delete z;
        return;
    }

    Node* l = z->left;
    Node* r = z->right;
    if (l == NIL && r == NIL) {
        if (p->left == z) p->left = NIL;
        else p->right = NIL;
        delete z;
    }
    else if (l == NIL || r == NIL) {
        Node* y = (l != NIL)? l : r;
        if (p->left == z) p->left = y;
        else p->right = y;
        y->parent = p;
        delete z;
    }
    else {
        Node* y = r;
        while (true) {
            if (y->left == NIL) break;
            y = y->left;
        }
        z->key = y->key;
        del(y);
    }
}

void del(int k)
{
    Node* z = find(k);
    if (z == NIL) return;
    del(z);
}

void inParse(Node* r)
{
    if (r == NIL) return;
    inParse(r->left);
    cout << " " << r->key;
    inParse(r->right);
}

void preParse(Node* r)
{
    if (r == NIL) return;
    cout << " " << r->key;
    preParse(r->left);
    preParse(r->right);
}

void print(Node* r) 
{
    inParse(r); cout << endl;
    preParse(r); cout << endl;
}

int main()
{
    int n;
    cin >> n;
    for (auto i = 0; i < n; i++) {
        char com[20];
        int k;
        cin >> com;
        if (com[0] == 'i') {
            cin >> k;
            insert(k);
        }
        else if (com[0] == 'p') {
            print(root);
        }
        else if (com[0] == 'f') {
            cin >> k;
            cout << ((find(k) != NIL)? "yes" : "no") << endl;
        }
        else if (com[0] == 'd') {
            cin >> k;
            del(k);
        }
    }

    return 0;
}
#endif

#ifdef ALDS1_7_D

struct Node {
    Node() : parent(-1), left(-1), right(-1) {};
    int parent;
    int left, right;
};

map<int, Node> nodes;

void postParse(int i)
{
    if (i == -1) return;
    postParse(nodes[i].left);
    postParse(nodes[i].right);
    static bool first = true;
    if (first) {
        first = false;
    } else {
        cout << " ";
    }
    cout << i;
}

void buildTree(
    int root, set<int>& remain, 
    map<int, int>& preorder, map<int, int>&inorder
)
{
    remain.erase(root);
    set<int> left, right;
    for (auto&& e : remain) {
        if (inorder[e] < inorder[root]) {
            left.insert(e);
        } else {
            right.insert(e);
        }
    }
    int ltop = -1;
    if (!left.empty()) {
        ltop = *min_element(left.begin(), left.end(), [&](int i, int j) {
            return preorder[i] < preorder[j];
        });
        buildTree(ltop, left, preorder, inorder);
    }
    int rtop = -1;
    if (!right.empty()) {
        rtop = *min_element(right.begin(), right.end(), [&](int i, int j) {
            return preorder[i] < preorder[j];
        });
        buildTree(rtop, right, preorder, inorder);
    }
    if (ltop != -1) {
        nodes[root].left = ltop;
        nodes[ltop].parent = root;
    }
    if (rtop != -1) {
        nodes[root].right = rtop;
        nodes[rtop].parent = root;
    }
}

int main()
{
    int n;
    cin >> n;
    vector<int> preOrd(n);
    vector<int> inOrd(n);
    for (auto i = 0; i < n; i++) {
        cin >> preOrd[i];
    }
    for (auto i = 0; i < n; i++) {
        cin >> inOrd[i];
    }

    map<int, int> preorder;  //up --> down position
    for (int i = 0; i < n; i++) {
        preorder[preOrd[i]] = i;
    }

    map<int, int> inorder;   //left --> right position
    for (int i = 0; i < n; i++) {
        inorder[inOrd[i]] = i;
    }

    int root = preOrd[0];
    set<int> remain(preOrd.begin(), preOrd.end());
    
    buildTree(root, remain, preorder, inorder);

    postParse(root); cout << endl;

    return 0;
}

#endif

#ifdef ALDS1_7_C
struct Node {
    Node() : parent(-1), left(-1), right(-1) {};
    int parent;
    int left, right;
};

map<int, Node> nodes;

void preParse(int i)
{
    if (i == -1) return;
    cout << " " << i;
    preParse(nodes[i].left);
    preParse(nodes[i].right);
}

void inParse(int i)
{
    if (i == -1) return;
    inParse(nodes[i].left);
    cout << " " << i;
    inParse(nodes[i].right);
}

void postParse(int i)
{
    if (i == -1) return;
    postParse(nodes[i].left);
    postParse(nodes[i].right);
    cout << " " << i;
}

int main()
{
    int n;
    cin >> n;
    for (auto i = 0; i < n; i++) {
        int id = -1, l = -1, r = -1;
        cin >> id >> l >> r;
        Node node;
        node.left = l;
        node.right = r;
        nodes[id] = node;
    }

    for (auto& it : nodes) {
        int id = it.first;
        const Node& node = it.second;
        if (node.left >= 0) {
            nodes[node.left].parent = id;
        }
        if (node.right >= 0) {
            nodes[node.right].parent = id;
        }
    }

    int root = -1;
    for (const auto& it : nodes) {
        const Node& node = it.second;
        if (node.parent == -1) {
            root = it.first;
            break;
        }
    }

    cout << "Preorder" << endl;
    preParse(root); cout << endl;

    cout << "Inorder" << endl;
    inParse(root); cout << endl;

    cout << "Postorder" << endl;
    postParse(root); cout << endl;

    return 0;
}
#endif

#ifdef ALDS1_6_B
struct Node {
    Node() : parent(-1), left(-1), right(-1) {};
    int parent;
    int left, right;
};

int GetType(const Node& N) {
    if (N.parent == -1) return 0; // root
    if (N.left == -1 && N.right == -1) return 1; // leaf
    else return 2;
}

string GetTypeStr(const Node& N) {
    int type = GetType(N);
    if (type == 0) return "root";
    else if (type == 1) return "leaf";
    else return "internal node";
}

int GetDegree(const Node& N) {
    int d = 0;
    if (N.left >= 0) d++;
    if (N.right >= 0) d++;
    return d;
}

map<int, Node> nodes;

int GetDepth(const Node& N) {
    if (GetType(N) == 0) return 0;
    else {
        return GetDepth(nodes[N.parent]) + 1;
    }
}

int GetSibling(int id) {
    const Node& n = nodes[id];
    if (n.parent >= 0) {
        const Node& p = nodes[n.parent];
        if ( p.left == id ) {
            return p.right;
        }
        else {
            return p.left;
        }
    }
    return -1;
}

int GetHeight(const Node& N) {
    if (N.left == -1 && N.right == -1) return 0;
    else {
        int lHeight = 0, rHeight = 0;
        if (N.left >= 0)
            lHeight = GetHeight(nodes[N.left]);
        if (N.right >= 0)
            rHeight = GetHeight(nodes[N.right]);
        return max(lHeight, rHeight) + 1;
    }
}


int main()
{
    int n;
    cin >> n;
    for (auto i = 0; i < n; i++) {
        int id = -1, l = -1, r = -1;
        cin >> id >> l >> r;
        Node node;
        node.left = l;
        node.right = r;
        nodes[id] = node;
    }

    for (auto& it : nodes) {
        int id = it.first;
        const Node& node = it.second;
        if (node.left >= 0) {
            nodes[node.left].parent = id;
        }
        if (node.right >= 0) {
            nodes[node.right].parent = id;
        }
    }

    for (const auto& it : nodes) {
        int id = it.first;
        const Node& node = it.second;
        cout << "node " << id << ": parent = " << node.parent;
        cout << ", sibling = " << GetSibling(id);
        cout << ", degree = " << GetDegree(node);
        cout << ", depth = " << GetDepth(node);
        cout << ", height = " << GetHeight(node);
        cout << ", " << GetTypeStr(node);
        cout << endl;
    }

    return 0;
}
#endif

#ifdef ALDS1_6_A
struct Node {
    Node() : parent(-1) {};
    int parent;
    vector<int> child;
};

map<int, Node> nodes;

int GetType(const Node& N) {
    if (N.parent == -1) return 0; // root
    if (N.child.empty()) return 1; // leaf
    else return 2;
}

string GetTypeStr(const Node& N) {
    int type = GetType(N);
    if (type == 0) return "root";
    else if (type == 1) return "leaf";
    else return "internal node";
}

int GetDepth(const Node& N) {
    if (GetType(N) == 0) return 0;
    else {
        return GetDepth(nodes[N.parent]) + 1;
    }
}

int main()
{
    int n;
    cin >> n;
    for (auto i = 0; i < n; i++) {
        int id = -1, k = 0;
        cin >> id >> k;
        Node node;
        for (auto j = 0; j < k; j++) {
            int ck = -1;
            cin >> ck;
            node.child.push_back(ck);
        }
        nodes[id] = node;
    }

    for (auto&& node : nodes) {
        for (auto&& cid : node.second.child ) {
            nodes[cid].parent = node.first;
        }
    }

    for (const auto& it : nodes) {
        int id = it.first;
        const Node& node = it.second;
        cout << "node " << id << ": parent = " << node.parent;
        cout << ", depth = " << GetDepth(node);
        cout << ", " << GetTypeStr(node);
        cout << ", [";
        for (auto k = 0; k < node.child.size(); k++) {
            if (k > 0) cout << ", ";
            cout << node.child[k];
        }
        cout << "]" << endl;
    }

    return 0;
}
#endif

#ifdef ALDS1_6_D
// [4, 3, 2, 7, 1, 6, 5] => 24
// [2, 1, 8, 10, 7, 9] => 49
int main()
{
    int n;
    cin >> n;
    vector<int> w(n);
    for (auto i = 0; i < n; i++) {
        cin >> w[i];
    }
    int wmin = *min_element(w.begin(), w.end());

    vector<int> ws = w;
    sort(ws.begin(), ws.end());

    map<int, int> order;
    for (auto i = 0; i < n; i++) {
        order[ws[i]] = i;
    }

    vector<bool> check(n, false);
    vector<set<int>>  cycle;
    for (auto i = 0; i < n; i++ ){
        set<int> c;
        int cur = i;
        while (true) {
            if (check[cur]) break;
            c.insert(w[cur]);
            check[cur] = true;
            cur = order[w[cur]];
        }
        if (c.size() > 1) {
            cycle.push_back(c);
        }
        /*
        for (auto&& e : c) {
            cout << " " << e; 
        }
        cout << endl;
        */
    }

    int cost = 0;
    for (auto&& c : cycle) {
        int wsum = accumulate(c.begin(), c.end(), 0);
        int cn = (int)c.size();
        int cmin = *c.begin();
        int cost1 = wsum + (cn - 2)*cmin; 
        int cost2 = wsum * 2;
        if (cmin > wmin) {
            cost2 = 2*(cmin + wmin) + wsum - cmin + wmin + (cn-2)*wmin;
        }
        cost += min(cost1, cost2);
        // cout << "cost1=" << cost1 << " " << "cost2=" << cost2 << endl;
    }

    cout << cost << endl;

    return 0;
}

#endif
#if 0

// ALDS1_6_C
template <class T>
int partition(T A[], int p, int r)
{
    auto x = A[r];
    auto i = p - 1;
    for (auto j = p; j < r; j++) {
        if (A[j] <= x) {
            i++;
            swap(A[i], A[j]);
        }
    }
    swap(A[i+1], A[r]);
    return i + 1;
}

template <class T>
void quickSort(T A[], int p, int r)
{
    if (p < r) {
        int q = partition<T>(A, p, r);
        quickSort(A, p, q-1);
        quickSort(A, q+1, r);
    }
}

template <class T>
void merge(T A[], int left, int mid, int right, T inf)
{
    int n1 = mid - left;
    int n2 = right - mid;

    vector<T> L(n1+1), R(n2+1);
    copy(A + left, A + mid, L.begin());
    copy(A + mid, A + right, R.begin());
    L[n1] = R[n2] = inf;

    int i = 0, j = 0;
    for (int k = left; k < right; k++ ) {
        if (L[i] <= R[j]) {
            A[k] = L[i];
            i++;
        } else {
            A[k] = R[j];
            j++;
        }
    }
}

template <class T>
void mergeSort(T A[], int left, int right, T inf)
{
    if (left + 1 < right) {
        int mid = (left + right)/2;
        mergeSort(A, left, mid, inf);
        mergeSort(A, mid, right, inf);
        merge(A, left, mid, right, inf);
    }
}

struct Card {
    char suit;
    int  number;
    bool operator<=(const Card& r) const { return !(this->number > r.number); }  
};

int main()
{
    int n;
    cin >> n;

    static Card A[100000];
    for (auto i = 0; i < n; i++) {
        cin >> A[i].suit >> A[i].number;
    }

    static Card B[100000];
    copy(A+0, A+n, B+0);

    quickSort<Card>(A, 0, n-1);

    Card INF = {'?', 1000000001};
    mergeSort<Card>(B, 0, n, INF);

    bool  stable = true;
    for (int i = 0; i < n; i++) {
        if (A[i].suit != B[i].suit) stable = false;
    }

    cout << (stable? "Stable" : "Not stable") << endl;

    for (auto i = 0; i < n; i++) {
        cout << A[i].suit << " " << A[i].number << endl;
    }
}

// ALDS1_6_B
int main()
{
    int n;
    cin >> n;

    static int A[100000];
    for (auto i = 0; i < n; i++) {
        cin >> A[i];
    }

    auto k = partition(A, 0, n-1);

    for (auto i = 0; i < n; i++) {
        if (i > 0) cout << " ";
        if (i == k) cout << "[";
        cout << A[i];
        if (i == k) cout << "]";
    }
    cout << endl;
}

// ALDS1_6_A
int main()
{
    const long int MAX = 2000000;
    long int n;  // 1 <= n <= MAX
    cin >> n;
    static array<int, MAX+1> A;  // 0 <= A[j] <= K, 0-origin
    for (auto i = 0; i < n; i++) {
        cin >> A[i];
    }
    const int K = 10000;

    static array<int, K+1> C; // counter
    for (auto i = 0; i <= K; i++ ) C[i] = 0;

    for (auto j = 0; j < n; j++) {
        C[A[j]]++;
    }

    // C[i] := number of A[j] <= i
    for (auto i = 1; i <= K; i++) {
        C[i] = C[i] + C[i-1];
    }

    static array<int, MAX+2> B; // sorted, 1-origin

    for (auto j = n-1; j >= 0; j--) {
        B[C[A[j]]] = A[j];
        C[A[j]]--;
    }

    for (auto i = 1; i <= n; i++) {
        if (i > 1) cout << " ";
        cout << B[i];
    }
    cout << endl;
}

// ALDS1_5_D
long int cnt = 0;
const int MAX = 500000;
long int L[MAX/2 +  2], R[MAX/2 + 2];

void merge( long int A[], int left, int mid, int right )
{
    int n1 = mid - left;
    int n2 = right - mid;

    copy(A + left, A + mid, L + 0);
    copy(A + mid, A + right, R + 0);
    L[n1] = R[n2] = LONG_MAX;

    int i = 0, j = 0;
    for (int k = left; k < right; k++ ) {
        if (L[i] <= R[j]) {
            A[k] = L[i];
            i++;
        } else {
            A[k] = R[j];
            j++;
            cnt += (n1 - i);
        }
    }
}

void mergeSort(long int A[], int left, int right)
{
    if (left + 1 < right) {
        int mid = (left + right)/2;
        mergeSort(A, left, mid);
        mergeSort(A, mid, right);
        merge(A, left, mid, right);
    }
}

int main()
{
    int n;
    cin >> n;
    vector<long int> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    mergeSort(&v[0], 0, n);

    cout << cnt << endl;
}

// ALDS1_5_C
const double PI=3.14159265358979323846;

struct Point {
    Point() { xy[0] = xy[1] = 0; }
    Point(double x, double y) { xy[0] = x; xy[1] = y; }
    double xy[2];
};

Point divide(const Point& p0, const Point& p1, double t)
{
    double new_xy[2];
    for (int k = 0; k < 2; k++) {
        new_xy[k] = (1. - t)*p0.xy[k] + t*p1.xy[k];
    }
    return Point(new_xy[0], new_xy[1]);
}

Point rotate(const Point& o, const Point& p, double angle)
{
    double RMat[2][2] = { { cos(angle), -sin(angle)}, {sin(angle), cos(angle)} };
    double vec[2] = {p.xy[0] - o.xy[0], p.xy[1] - o.xy[1]};
    double vec2[2];
    for (int k = 0; k < 2; k++) {
        vec2[k] = RMat[k][0]*vec[0] + RMat[k][1]*vec[1];
    }
    return Point(o.xy[0] + vec2[0], o.xy[1] + vec2[1]);
}

void kock_divide( const Point& p0, const Point& p1, Point mid[3] )
{
    mid[0] = divide(p0, p1, 1/3.0);
    mid[2] = divide(p0, p1, 2/3.0);
    mid[1] = rotate(mid[0], mid[2], PI/3.0);
}

void kock(const Point& p0, const Point& p1, int n, vector<Point>* vv)
{
    vv->clear();
    if (n == 0) {
        vv->push_back(p0);
        vv->push_back(p1);
        return;
    }
    Point mid[3];
    kock_divide(p0, p1, mid);

    vector<Point> v[4];
    kock(p0, mid[0], n-1, &v[0]);
    kock(mid[0], mid[1], n-1, &v[1]);
    kock(mid[1], mid[2], n-1, &v[2]);
    kock(mid[2], p1, n-1, &v[3]);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < v[i].size()-1; j++) {
            vv->push_back( v[i][j] );
        }
    }
    vv->push_back(p1);
}

int main()
{
    int n;
    cin >> n;

    vector<Point> v;
    kock(Point(0,0), Point(100,0), n, &v);

    cout << fixed << setprecision(8); 
    for (auto i = 0; i < v.size(); i++) {
        cout << v[i].xy[0] << " " << v[i].xy[1] << endl;
    }
}

// ALDS1_5_B
int cnt = 0;
const int MAX = 500000;
long int L[MAX/2 +  2], R[MAX/2 + 2];

void merge( long int A[], int left, int mid, int right )
{
    int n1 = mid - left;
    int n2 = right - mid;

    copy(A + left, A + mid, L + 0);
    copy(A + mid, A + right, R + 0);
    L[n1] = R[n2] = LONG_MAX;

    int i = 0, j = 0;
    for (int k = left; k < right; k++ ) {
        cnt++;
        if (L[i] <= R[j]) {
            A[k] = L[i];
            i++;
        } else {
            A[k] = R[j];
            j++;
        }
    }
}

void mergeSort(long int A[], int left, int right)
{
    if (left + 1 < right) {
        int mid = (left + right)/2;
        mergeSort(A, left, mid);
        mergeSort(A, mid, right);
        merge(A, left, mid, right);
    }
}

int main()
{
    int n;
    cin >> n;
    vector<long int> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    mergeSort(&v[0], 0, n);

    for (int i = 0; i < n; i++) {
        if (i > 0) cout << " ";
        cout << v[i];
    }
    cout << endl;
    cout << cnt << endl;
}

// ALDS1_5_A
bool sumUpTo(const int v[], int n, int m)
{
    if (n == 1) {
        if (v[0] == m) return true;
        if (m == 0) return true;
        else return false;
    }
    if (sumUpTo(v+1, n-1, m-v[0])) {
        return true;
    }
    if (sumUpTo(v+1, n-1, m)) {
        return true;
    }
    return false;
}

int main()
{
    int n;
    cin >> n;
    vector<int> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int m;
        cin >> m;

        if (sumUpTo(&v[0], n, m)) {
            cout << "yes" << endl;
        } else {
            cout << "no" << endl;
        }
    }
}

// ALDS1_4_D
int calcTrackNumber(const vector<int>& w, int p)
{
    int count = 1;
    int rest = p;
    for (auto i = 0; i < w.size(); i++) {
        if (w[i] <= rest) {
            rest -= w[i];
        }
        else {
            count++;
            rest = p - w[i];
        }
    }
    return count;
}

int main()
{
    int n, k;
    cin >> n >> k;
    vector<int> w(n);
    int wSum = 0;
    int wMax = 0;
    for (int i = 0; i < n; i++) {
        cin >> w[i];
        wSum += w[i];
        wMax = max(wMax, w[i]);
    }

    int s = wMax - 1;  // Ng
    int t = wSum;      // Ok

    while (t - s > 1) { 
        int m = s + (t - s)/2;
        int num = calcTrackNumber(w, m);
        if (num <= k) { // m is Ok
            t = m;
        }
        else {  // m is Ng
            s = m;
        }
    }
    cout << t << endl;

    return 0;
}

// ALDS1_4_C

int charInt(char c)
{
    if (c == 'A') return 1;
    if (c == 'C') return 2;
    if (c == 'G') return 3;
    if (c == 'T') return 4;
    return 0;
}

int strKey(const char str[])
{
    int s = 0, p = 1;
    for (auto i = 0; i < strlen(str); i++, p*= 5 ){
        s += p * charInt(str[i]);
    }
    return s;
}

const int M = 1046527;  // Carol prime

int Hash1(int k) 
{
    return k % M;
}

int Hash2(int k) 
{
    return 1 + k % (M -1);
}

int Hash(int k, int i)
{
    return (Hash1(k) + i*Hash2(k)) % M;
}

array<char[13], M> hashTable;

const int NIL = -1;

void init()
{
    for (int i = 0; i < M; i++) hashTable[i][0] = '\0';
}

int find(const char str[])
{
    int k = strKey(str);
    for (int i = 0;; i++) {
        int h = Hash(k, i);
        if (strcmp(hashTable[h], str) == 0) {
            return 1;
        }
        else if (hashTable[h][0] == '\0') {
            return 0;
        }
    }
    return 0;
}

int insert(const char str[])
{
    int k = strKey(str);
    for (int i = 0;; i++) {
        int h = Hash(k, i);
        if (strcmp(hashTable[h], str) == 0) {
            return 1;
        }
        else if (hashTable[h][0] == '\0') {
            strcpy(hashTable[h], str);
            return 0;
        }
    }
    return 0;
}

int main()
{
    init();

    int n;
    cin >> n; 
    for (int i = 0; i < n; i++) {
        char com[20];
        char s[13];
        cin >> com >> s;

        if (com[0] == 'i') {
            insert(s);
        }
        else if (com[0] == 'f') {
            cout << (find(s)? "yes" : "no") << endl;
        }
    }

    return 0;
}

// ALDS1_4_B
int main()
{
    static array<int, 100001> a;
    int n;
    cin >> n; 
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    sort(a.begin(), a.begin()+n);

    int count = 0;
    int m;
    cin >> m;
    for (int i = 0; i < m; i++) {
        int x;
        cin >> x;
        int s = 0;
        int t = n-1;
        while (t - s >= 0) {
            int m = s + (t - s)/2;
            if (a[m] == x) {
                count++;
                break;
            }
            else if (a[m] < x) {
                s = m + 1;
            }
            else if (a[m] > x) {
                t = m - 1;
            }
        }
    }

    cout << count << endl;

    return 0;
}

// ALDS1_4_A
int main()
{
    static array<int, 10001> a;
    int n;
    cin >> n; 
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    int count = 0;
    int m;
    cin >> m;
    for (int i = 0; i < m; i++) {
        int x;
        cin >> x;
        if ( find(a.begin(), a.begin()+n, x) != a.begin()+n) count++;
    }

    cout << count << endl;

    return 0;
}

// ALDS1_3_D
int main()
{
    string str;
    getline(cin, str);

    stack<int> s;
    vector<pair<int, int>> s2;

    int sum = 0;
    int n = (int)str.size();
    for (auto i = 0; i < n; i++ ) {
        if (str[i] == '\\') {
            s.push(i);
        }
        else if (str[i] == '/' && !s.empty()) {
            int j = s.top();
            int area = i - j;
            sum += area;
            s.pop();
            while (!s2.empty() && j < s2.back().first) {
                area += s2.back().second;
                s2.pop_back();
            }
            s2.push_back(make_pair(j, area));
        }
    }

    cout << sum << endl;
    cout << s2.size();
    for (auto i = 0; i < s2.size(); ++i) {
        cout << " " << s2[i].second;
    }
    cout << endl;

    return 0;
}

// ALDS1_3_C
struct Node {
    Node(int x) : key(x), prev(nullptr), next(nullptr) {}
    int key;
    Node* prev;
    Node* next;
};

class List
{
public:
    List() {
        nil = new Node(0);
        nil->prev = nil;
        nil->next = nil;
    }
    ~List() {
        free();
    }
    void insert(int x) {
        Node* m = new Node(x);
        nil->next->prev = m;
        m->next = nil->next;
        m->prev = nil;
        nil->next = m;
    }
    Node* find(int x) {
        Node* cur = nil->next;
        while ( cur != nil && cur->key != x) {
            cur = cur->next;
        }
        return cur;
    }
    void remove(Node* m) {
        if (m == nullptr) return;
        if (m == nil) return;
        m->prev->next = m->next;
        m->next->prev = m->prev;
        delete m;
    }
    Node* front() {
        if (nil->next != nil) {
            return nil->next;
        } else {
            return nullptr;
        }
    }
    Node* back() {
        if (nil->prev != nil) {
            return nil->prev;
        } else {
            return nullptr;
        }
    }
    void free() {
        Node* cur = nil->next;
        while ( cur != nil ) {
            Node* next = cur->next;
            delete cur;
            cur = next;
        }
        delete nil;
    }

private:
    Node* nil;
};

int main()
{
    List L;

    int n;
    scanf("%d", &n);

    for (int i = 0; i < n; i++ ) {
        char com[20];
        int  x;
        scanf("%s", com);
        if (com[0] == 'i') {
            scanf("%d", &x);
            // insert x
            L.insert(x);
        }
        else if (com[0] == 'd') {
            if (strlen(com) > 6 && com[6] == 'F') {
                // deleteFirst
                L.remove(L.front());
            }
            else if (strlen(com) > 6 && com[6] == 'L') {
                // deleteLast
                L.remove(L.back());
            }
            else {
                scanf("%d", &x);
                // delete x
                L.remove(L.find(x));
            }
        }
    }

    int cnt = 0;
    Node* cur = L.front();
    while (cur) {
        if (cnt > 0) printf(" ");
        printf("%d", cur->key);
        cnt++;
        if (cur == L.back()) break;
        cur = cur->next;
    }
    printf("\n");


    return 0;
}

// ALDS1_3_B
template<class T, int MAX = 100>
class Queue
{
public:
    Queue<T, MAX>() : head(0), tail(0) {}
    T pop_front() {
        if (empty()) throw "No element";
        T x = Q[head++];
        if (head == MAX) head = 0;
        return x;
    } 
    void push_back(T x) {
        if (full()) throw "Full";
        Q[tail++] = x;
        if (tail == MAX) tail = 0;
    }
    bool empty() const {
        return (tail - head) % MAX == 0;
    }
    bool full() const {
        return (tail - head) % MAX == MAX-1;
    }
private:
    int head;
    int tail;
    array<T, MAX> Q;
};

struct Proc {
    char name[11];
    int  t;
};

int main()
{
    static Queue<Proc, 100001> Q;

    int n, q;
    cin >> n >> q;

    for (int i = 0; i < n; i++ ) {
        Proc p;
        cin >> p.name >> p.t; 
        Q.push_back(p);
    }

    int t = 0;
    while (!Q.empty()) {
        auto p = Q.pop_front();
        if (p.t > q) {
            p.t -= q;
            Q.push_back(p);
            t += q;
        } else {
            t += p.t;
            cout << p.name << " " << t << endl;
        }
    }

    return 0;
}

// ALDS1_3_A
template<class T>
class Stack
{
public:
    T pop() {
        if (n < 1) throw "No element";
        return S[--n];
    } 
    void push(T x) {
        if (n >= 99) throw "Full";
        S[n++] = x;
    }
private:
    int n = 0;
    array<T, 100> S;
};

int main()
{
    Stack<int> S;
    string str;
    while (cin >> str) {
        if (str == "+" || str == "*" || str == "-") {
            int b = S.pop();
            int a = S.pop();
            int res;
            if (str == "+" ) res = a + b;
            if (str == "*" ) res = a * b;
            if (str == "-" ) res = a - b;
            S.push(res);
        }
        else {
            S.push(atoi(str.c_str()));
        }
     }
    cout << S.pop() << endl;

    return 0;
}

template<class T, size_t N>
void print(const array<T, N>& a, int n)
{
    for (int k = 0; k < n; k++) {
        if (k > 0) cout << " ";
        cout << a[k];
    }
    cout << endl;
}

// ALDS1_2_D
void insertionSort(int a[], int n, int g, long long& cnt)
{
    for (int i = g; i < n; i++) {
        int v = a[i];
        int j = i - g;
        while (j >= 0 && a[j] > v) {
            a[j + g] = a[j];
            j -= g;
            cnt++;
        }
        a[j + g] = v;
    }
}

int main()
{
    array<int, 1000000> a;

    int n;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    vector<int> G;  // 1, 4, 13, 40, ... (g(n+1) = 3g(n) + 1
    int g = 0;
    for (int g = 1; ; ) {
        if (g > n) break;
        G.push_back(g);
        g = 3*g + 1;
    }

    long long cnt = 0;
    for (auto it = G.rbegin(); it != G.rend(); ++it) {
        insertionSort(&a[0], n, *it, cnt);
    }

    cout << G.size() << endl;
    for (int i = (int)G.size()-1; i >= 0; i--) {
        cout << G[i];
        if (i > 0) cout << " ";
    }
    cout << endl;

    cout << cnt << endl;

    for (int i = 0; i < n; i++) {
        cout << a[i] << endl;
    }

    return 0;
}

// ALDS1_2_C
struct Card {
    char suit; 
    int  value;
};

template<size_t N>
void print(const array<Card,N>& a, int n)
{
    for (int k = 0; k < n; k++) {
        if (k > 0) cout << " ";
        cout << a[k].suit << a[k].value;
    }
    cout << endl;
}

void bubbleSort(Card a[], int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = n-1;  j > i; j--) {
            if (a[j].value < a[j-1].value )
                swap(a[j], a[j-1]);
        }
    }
}

void selectionSort(Card a[], int n)
{
    for (int i = 0; i < n; i++) {
        int iMin = i;
        for (int j = i; j < n; j++) {
            if(a[j].value < a[iMin].value) 
                iMin = j;
        }
        swap(a[i], a[iMin]);
    }
}

int main()
{
    array<Card, 100> a, b;

    int n;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> a[i].suit >> a[i].value;
    }
    b = a;

    bubbleSort(&a[0], n);
    print(a, n);
    cout << "Stable" << endl;

    selectionSort(&b[0], n);
    print(b, n);

    bool  stable = true;
    for (int i = 0; i < n; i++) {
        if (a[i].suit != b[i].suit) stable = false;
    }

    cout << (stable? "Stable" : "Not stable") << endl;

    return 0;
}

// ALDS1_2_B
int main()
{
    array<int, 100> a;

    int n;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    int sw = 0;
    for (int i = 0; i < n; i++) {
        int iMin = i;
        for (int j = i; j < n; j++ ) {
            if (a[j] < a[iMin]) {
                iMin = j;
            }
        }
        if (i != iMin) {
            sw++;
            swap(a[i], a[iMin]);
        }
    }

    print(a, n);
    cout << sw << endl;

    return 0;
}

// ALDS1_2_A
int main()
{
    array<int, 100> a;

    int n;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    int sw = 0;
    bool flag = true;
    int i = 0;
    while (flag) {
        flag = false;
        for (int j = n-1; j > i; j--) {
            if (a[j] < a[j-1]) {
                swap(a[j], a[j-1]);
                flag = true;
                sw++;
            }
        }
        i++;
    }

    print(a, n);
    cout << sw << endl;

    return 0;
}

// ALDS1_1_D
int main()
{
    int n;
    cin >> n;
    int min;
    cin >> min;
    int r_max = INT_MIN;
    for (int i = 1; i < n; i++) {
        int x;
        cin >> x;

        if (x - min > r_max) {
            r_max = x - min;
        }
        if (x < min) min = x;
    }
    cout << r_max << endl;

    return 0;
}

// ALDS1_1_C
bool isPrime(int x) {
    for (int i = 2; i*i <= x; i++) {
        if (x % i == 0)
            return false;
    }
    return true;
}

int main()
{
    int n;
    cin >> n;
    int count = 0;
    for (int i = 0; i < n; i++) {
        int x;
        cin >> x;

        if (isPrime(x)) {
            count++;
        }
    }
    cout << count << endl;

    return 0;
}

// ALDS1_1_B
int main()
{
    long long x, y;
    cin >> x >> y;

    if (x > y) swap(x, y);
    while (y > 0) {
        if (x > y) swap(x, y);
        y = y % x;
    }
    cout << x << endl;

    return 0;
}

// ALDS1_1_A
int main()
{
    array<int, 100> a;

    int n;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }
    print(a, n);
    for (int i = 1; i < n; i++ ) {
        int v = a[i];
        int j = i-1;
        while (j >= 0 && a[j] > v) {
            a[j+1] = a[j];
            j--;
        }
        a[j+1] = v;

        print(a, n);
    }

    return 0;
}

#endif
