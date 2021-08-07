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

int main()
{
    return 0;
}

#ifdef GRL_6_B

const int INF = 1 << 21;

struct graph {
    explicit graph(int n) { g = vector<vector<edge>>(n); }

    void add_edge(int from, int to, int c, int d) {
        g[from].push_back(edge(to, c, d, (int)g[to].size()));  
        g[to].push_back(edge(from, 0, -d, (int)g[from].size()-1));
    }

    int min_cost_flow(int s, int t, int f)
    {
        int res = 0;
        vector<int> d(g.size());
        vector<pair<int, int>> prev(g.size());
        while (f > 0) {
            bool valid = bellman_ford(s, d, prev);
            if (!valid) {
                return -1;
            }
            if (d[t] == INF) {
                return -1;
            }

            int minf = f;
            for (int v = t; v != s; v = prev[v].first) {
                minf = min(minf, g[prev[v].first][prev[v].second].capacity);
            }
            f -= minf;
            res += minf * d[t];
            for (int v = t; v != s; v = prev[v].first) {
                edge& e = g[prev[v].first][prev[v].second];
                e.capacity -= minf;
                g[v][e.rev].capacity += minf;
            }
        }
        return res;
    }

    bool bellman_ford(int s, vector<int>& d, vector<pair<int, int>>& prev)
    {
        fill(d.begin(), d.end(), INF);
        d[s] = 0;

        int cnt = 0;
        bool updated = true;
        while (updated) {
            updated = false;
            for ( int v = 0; v < (int)g.size(); v++ ) {
                if (d[v] == INF) continue;
                for ( int i = 0; i < (int)g[v].size(); i++ ) {
                    edge& e = g[v][i];
                    if (e.capacity > 0) {
                        if (d[e.to] > d[v] + e.dst) {
                            d[e.to] = d[v] + e.dst;
                            prev[e.to] = make_pair(v, i);
                            updated = true;
                        }
                    }
                }
            }
            cnt++;
            if (cnt == (int)g.size() && cnt > 2) {
                return false;
            }
        }
        return true;
    }

private:
    struct edge {
        int to;
        int capacity;  // rsidual capacity
        int dst;       // cost
        int rev;    // reverse-edge id in g[to]
        edge(int t, int c, int d, int r) : to(t), capacity(c), dst(d), rev(r) {}
    };

    vector<vector<edge>> g;
};

int main()
{
    int n, e, f;
    cin >> n >> e >> f;

    graph G(n);

    for (int i = 0; i < e; i++) {
        int from, to, cap, dst;
        cin >> from >> to >> cap >> dst;
        G.add_edge(from, to, cap, dst);
    }

    int minc = G.min_cost_flow(0, n-1, f);

    cout << minc << endl;

    return 0;
}

#endif

#ifdef GRL_7_A

const int INF = 1 << 21;

struct graph {
    struct edge {
        int to;
        int capacity;  // rsidual capacity
        int rev;    // reverse-edge id in g[to]
        edge(int t, int c, int r) : to(t), capacity(c), rev(r) {}
    };

    vector<vector<edge>> g;

    explicit graph(int n) { g = vector<vector<edge>>(n); }

    void add_edge(int from, int to, int c) {
        g[from].push_back(edge(to, c, (int)g[to].size()));  
        g[to].push_back(edge(from, 0, (int)g[from].size()-1));
    }

    // Ford-Fulkerson
    int max_flow(int s, int t)
    {
        int flow = 0;
        while (true) {
            vector<bool>  visited(g.size(), false);
            int f = dfs(s, t, visited, INF);
            if (f == 0) {
                return flow;
            }
            flow += f;
        }
    }

private:
    int dfs(int v, int t, vector<bool>& visited, int f)
    {
        if (v == t) return f;

        visited[v] = true;

        for (auto& e : g[v]) {
            if (!visited[e.to] && e.capacity > 0) {
                int fmin = min(f, e.capacity);
                int d = dfs(e.to, t, visited, fmin);
                if (d > 0) {
                    e.capacity -= d;
                    g[e.to][e.rev].capacity += d;
                    return d;
                }
            }
        }

        return 0;
    }
};

int main()
{
    int n, m, e;
    cin >> n >> m >> e;

    graph G(n + m + 2);

    for (int i = 0; i < e; i++) {
        int from, to;
        cin >> from >> to;
        G.add_edge(from, n + to, 1);
    }
    int s = n + m;
    int t = n + m + 1;

    for (int i = 0; i < n; i++) {
        G.add_edge(s, i, 1);
    }

    for (int i = n; i < n + m; i++) {
        G.add_edge(i, t, 1);
    }

    int flow = G.max_flow(s, t);

    cout << flow << endl;

    return 0;
}

#endif

#ifdef GRL_6_A

const int INF = 1 << 21;

struct edge {
    int to;
    int capacity;  // rsidual capacity
    int rev;    // reverse-edge id in g[to]
    edge(int t, int c, int r) : to(t), capacity(c), rev(r) {}
};

int dfs(int v, int t, vector<vector<edge>>& g, vector<bool>& visited, int f)
{
    if (v == t) return f;

    visited[v] = true;

    for (auto& e : g[v]) {
        if (!visited[e.to] && e.capacity > 0) {
            int fmin = min(f, e.capacity);
            int d = dfs(e.to, t, g, visited, fmin);
            if (d > 0) {
                e.capacity -= d;
                g[e.to][e.rev].capacity += d;
                return d;
            }
        }
    }

    return 0;
}

// Ford-Fulkerson
int max_flow(int s, int t, vector<vector<edge>>& g)
{
    int flow = 0;
    while (true) {
        vector<bool>  visited(g.size(), false);
        int f = dfs(s, t, g, visited, INF);
        if (f == 0) {
            return flow;
        }
        flow += f;
    }
}

int main()
{
    int n, e;
    cin >> n >> e;

    vector<vector<edge>> g(n);

    for (int i = 0; i < e; i++) {
        int from, to, c;
        cin >> from >> to >> c;
        g[from].push_back(edge(to, c, (int)g[to].size()));  
        g[to].push_back(edge(from, 0, (int)g[from].size()-1));
    }

    int flow = max_flow(0, n-1, g);

    cout << flow << endl;

    return 0;
}

#endif

#ifdef GRL_5_A

#define INF (1 << 21)

struct edge {
    int to;
    int cost;
    edge( int t, int w ) : to(t), cost(w) {}
};

void bfs(int s, const vector<vector<edge>>& g, vector<int>& dist)
{
    queue<int> Q;
    Q.push(s);
    dist[s] = 0;

    while (!Q.empty()) {
        int u = Q.front(); Q.pop();
        for (auto& e : g[u]) {
            if (dist[e.to] == INF) {
                dist[e.to] = dist[u] + e.cost;
                Q.push(e.to);
            }
        }
    }
}

int calc_diameter(const vector<vector<edge>>& g)
{
    int n = (int)g.size();

    vector<int> dist(n, INF);
    bfs(0, g, dist);

    int far = -1;
    int maxdist = 0;
    for (int i = 0; i < n; i++) {
        if (dist[i] != INF && maxdist < dist[i]) {
            maxdist = dist[i];
            far = i;
        }
    }

    fill(dist.begin(), dist.end(), INF);
    bfs(far, g, dist);

    maxdist = 0;
    for (int i = 0; i < n; i++) {
        if (dist[i] != INF && maxdist < dist[i]) {
            maxdist = dist[i];
            far = i;
        }
    }

    return maxdist;
}

int main()
{
    int n;
    cin >> n;

    vector<vector<edge>> g(n);

    for (int i = 0; i < n-1; i++) {
        int s, t, w;
        cin >> s >> t >> w;
        g[s].push_back(edge(t, w));
        g[t].push_back(edge(s, w));
    }

    int d = calc_diameter(g);

    cout << d << endl;

    return 0;
}

#endif

#ifdef GRL_4_B

void dfs(int s, const vector<vector<int>>& g, vector<bool>& visited, vector<int>& out)
{
    visited[s] = true;
    for (auto& v : g[s]) {
        if (!visited[v]) {
            dfs(v, g, visited, out);
        }
    }
    // posorder
    out.push_back(s);
}

void topo_sort_dfs(const vector<vector<int>>& g, vector<int>& out)
{
    int n = (int)g.size();
    vector<bool> visited(n, false);

    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            dfs(i, g, visited, out);
        }
    }

    reverse(out.begin(), out.end());
}

int main()
{
    int n, e;
    cin >> n >> e;

    vector<vector<int>> g(n);

    for (int i = 0; i < e; i++) {
        int s, t;
        cin >> s >> t;
        g[s].push_back(t);
    }

    vector<int> v;
    topo_sort_dfs(g, v);

    for (auto& u : v) {
        cout << u << endl;
    }

    return 0;
}

#endif

#ifdef GRL_4_A

void bfs(int s, const vector<vector<int>>& g, vector<bool>& visited, vector<int>& indeg, vector<int>& out)
{
    queue<int> Q;
    Q.push(s);

    while(!Q.empty()) {
        int u = Q.front(); Q.pop();
        visited[u] = true;
        out.push_back(u);
        for (auto & v : g[u]) {
            if (!visited[v]) {
                indeg[v]--;
                if (indeg[v] == 0) {
                    Q.push(v);
                }
            }
        }
    }
}

bool topo_sort_bfs(const vector<vector<int>>& g, vector<int>& out)
{
    int n = (int)g.size();
    vector<bool> visited(n, false);
    vector<int> indeg(n, 0);

    for (int i = 0; i < n; i++) {
        for (auto& v : g[i]) {
            indeg[v]++;
        }
    }

    for (int i = 0; i < n; i++) {
        if (!visited[i] && indeg[i] == 0) {
            bfs(i, g, visited, indeg, out);
        }
    }

    if ((int)out.size() == n) {
        return true;
    }
    else {
        return false;
    }
}

int main()
{
    int n, e;
    cin >> n >> e;

    vector<vector<int>> g(n);

    for (int i = 0; i < e; i++) {
        int s, t;
        cin >> s >> t;
        g[s].push_back(t);
    }

    vector<int> v;
    bool isDAG = topo_sort_bfs(g, v);

    cout << (isDAG? 0 : 1) << endl;

    return 0;
}

#endif

#ifdef GRL_3_C

void dfs(int r, int p, const vector<vector<int>>& g, vector<bool>& visited, vector<int>& postorder)
{
    visited[r] = true;
    for ( auto & v : g[r]) {
        if (v == p) {  
            continue;
        }
        if (visited[v]) { 
            continue;
        }

        dfs(v, r, g, visited, postorder);
    }
    postorder.push_back(r);
}

void r_dfs(int r, int p, const vector<vector<int>>& g, vector<bool>& visited, set<int>& group)
{
    visited[r] = true;
    group.insert(r);
    for ( auto & v : g[r]) {
        if (v == p) {  
            continue;
        }
        if (visited[v]) { 
            continue;
        }

        r_dfs(v, r, g, visited, group);
    }
}

vector<int> strongly_connected_group(const vector<vector<int>> g, const vector<vector<int>> rg)
{
    size_t v = g.size();

    vector<bool> visited(v, false);
    vector<int> postorder;

    for (int i = 0; i < v; i++) {
        int next = i;
        if (!visited[next]) {
            dfs(next, -1, g, visited, postorder);
        }
    }

    reverse(postorder.begin(), postorder.end());

    fill(visited.begin(), visited.end(), false);

    vector<int> sc_group(v, -1);
    int gr_id = 0;
    for (int i = 0; i < v; i++) {
        set<int> group;
        int next = postorder[i];
        if (!visited[next]) {
            r_dfs(next, -1, rg, visited, group);
            if (group.size() > 1) {
                for (auto& u : group) {
                    sc_group[u] = gr_id;
                }
                gr_id++;
            }
        }
    }

    return sc_group;
}

int main()
{
    int v, e;
    cin >> v >> e;

    vector<vector<int>> g(v);
    vector<vector<int>> rg(v);

    for (int i = 0; i < e; i++) {
        int s, t;
        cin >> s >> t;
        g[s].push_back(t);
        rg[t].push_back(s);
    }

    vector<int> sc_group = strongly_connected_group(g, rg);

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int u, v;
        cin >> u >> v;
        if (sc_group[u] >= 0 && sc_group[u] == sc_group[v] ) {
            cout << 1 << endl;
        } else {
            cout << 0 << endl;
        }
    }

    return 0;
}

#endif

#ifdef GRL_3_A

void dfs(int r, const vector<vector<int>>& g, vector<int>& parent, vector<int>& preorder, vector<int>& lowest)
{
    static int time = 0;
    preorder[r] = time++;
    lowest[r] = min(lowest[r], preorder[r]);
    for ( auto & v : g[r]) {
        if (v == parent[r]) {  
            continue;
        }
        if (preorder[v] >= 0) { // backedge (T‚É‘®‚³‚È‚¢G‚ÌƒGƒbƒW)
            lowest[r] = min(lowest[r], preorder[v]);
            continue;
        }

        parent[v] = r;
        dfs(v, g, parent, preorder, lowest);
        lowest[r] = min(lowest[r], lowest[v]);
    }
}

void articulation_point(const vector<vector<int>>& g, set<int>& atp)
{
    size_t v = g.size();
    vector<int> parent(v, -1);
    vector<int> preorder(v, -1);
    vector<int> lowest(v, (int)v);

    dfs(0, g, parent, preorder, lowest);

    int n0 = 0;
    for (int i = 1; i < v; i++) {
        int p = parent[i];
        if (p == 0) n0++;
        else if (preorder[p] <= lowest[i]) {
            atp.insert(p);
        }
    }
    if (n0 > 1) atp.insert(0);
}

int main()
{
    int v, e;
    cin >> v >> e;

    vector<vector<int>> g(v);

    for (int i = 0; i < e; i++) {
        int s, t;
        cin >> s >> t;
        g[s].push_back(t);
        g[t].push_back(s);
    }

    set<int> atp;
    articulation_point(g, atp);

    for (auto& p : atp) {
        cout << p << endl;
    }

    return 0;
}

#endif

#ifdef GRL_2_A

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

struct edge {
    int from;
    int to;
    int cost;
};

void kruskal(int v, vector<edge>& edges, vector<edge>& mst)
{
    sort(edges.begin(), edges.end(), 
        [](const auto& a, const auto& b) { return a.cost < b.cost; }
    );

    mst.clear();

    DisjointSet ds(v);
    for ( auto& e : edges ) {
        if (ds.findSet(e.from) != ds.findSet(e.to)) {
            ds.unite(e.from, e.to);
            mst.push_back(e);
        }
    }
}

int main()
{
    int v, e;
    cin >> v >> e;

    vector<edge> edges;

    for (int i = 0; i < e; i++) {
        int s, t, w;
        cin >> s >> t >> w;
        edge E = {s, t, w};
        edges.push_back(E);
    }

    vector<edge> mst;
    kruskal(v, edges, mst);

    long long sum = accumulate(mst.begin(), mst.end(), 
        0LL, 
        [](long long sum, auto& e){ return sum + (long long)e.cost;}
    );

    cout << sum << endl;

    return 0;
}

#endif

#ifdef GRL_1_C

using ll = long long;

static const int MAX = 100;
static const ll INF = LLONG_MAX;

ll dp[MAX][MAX];   // shortest peth length from [i] to [j]

int main()
{
    int n, e;
    cin >> n >> e;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) dp[i][j] = 0;
            else dp[i][j] = INF;
        }
    }

    for (int i = 0; i < e; i++) {
        int s, t, d;
        cin >> s >> t >> d;
        dp[s][t] = d;
    }

    // warshall-floyd
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            if (dp[i][k] == INF) continue;
            for (int j = 0; j < n; j++) {
                if (dp[k][j] == INF) continue;
                dp[i][j] = min(dp[i][j], dp[i][k] + dp[k][j]);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        if (dp[i][i] < 0) {
            cout << "NEGATIVE CYCLE" << endl;
            return 0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j) cout << " ";
            if (dp[i][j] < INF) {
                cout << dp[i][j];
            } else {
                cout << "INF";
            }
        }
        cout << endl;
    }

    return 0;
}

#endif

#ifdef GRL_1_B

using ll = long long;
const ll INF = LLONG_MAX;

struct edge {
    int from;
    int to;
    ll  cost;
};

// ‚Ç‚±‚©‚É•‰‚Ì•Â˜H‚ª‚ ‚é‚©
bool find_nagative_loop(int n, vector<edge>& g)
{
    vector<ll> d(n, 0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < (int)g.size(); j++) {
            edge& e = g[j];
            if (d[e.to] > d[e.from] + e.cost) {
                d[e.to] = d[e.from] + e.cost;
                if (i == n - 1) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool bellman_ford(int s, int n, vector<edge>& g, vector<ll>& d)
{
    fill(d.begin(), d.end(), INF);
    d[s] = 0;

    int cnt = 0;
    while (true) {
        bool updated = false;
        for (int j = 0; j < (int)g.size(); j++) {
            edge& e = g[j];
            if (d[e.from] <INF && d[e.to] > d[e.from] + e.cost) {
                d[e.to] = d[e.from] + e.cost;
                updated = true;
            }
        }
        if (!updated) {
            break;
        }
        cnt++;
        if (cnt == n) {
            return false;
        }
    }
    return true;
}

int main()
{
    int n, e, r;
    cin >> n >> e >> r;

    vector<edge> g;
    for (int i = 0; i < e; i++) {
        int s, t, d;
        cin >> s >> t >> d;
        edge E = {s, t, d};
        g.push_back(E);
    }

    vector<ll> d(n);
    bool b = bellman_ford(r, n, g, d);
    if (b == false) {
        cout << "NEGATIVE CYCLE" << endl;
        return 0;
    }

    for ( int i = 0; i < n; i++ ) {
        if (d[i] < INF) {
            cout << d[i] << endl;
        } else {
            cout << "INF" << endl;
        }
    }

    return 0;
}

#endif

#ifdef GRL_1_A

using ll = long long;
const ll INF = LLONG_MAX;

struct edge {
    int to;
    ll cost;
};
typedef pair<ll, int> P; // first:dist, second:vid

vector<vector<edge>> g;

vector<ll>  dijkstra(int s, int n)
{
    vector<ll> d(n, INF);

    d[s] = 0;
    priority_queue<P, vector<P>, greater<P>> PQ; 
    PQ.push(make_pair(0, s));

    while(!PQ.empty()) {
        auto cur = PQ.top(); PQ.pop();
        int u = cur.second;
        if (d[u] < cur.first) continue;
        for (auto& e : g[u]) {
            if (d[u] + e.cost < d[e.to]) {
                d[e.to] = d[u] + e.cost;
                PQ.push(make_pair(d[e.to], e.to));
            }
        }    
    }

    return d;
}

int main()
{
    int n, e, r;
    cin >> n >> e >> r;
    g.resize(n);

    for (int i = 0; i < e; i++) {
        int s, t, d;
        cin >> s >> t >> d;
        edge E = {t, d};
        g[s].push_back(E);
    }

    vector<ll> d = dijkstra(r, n);

    for ( int i = 0; i < n; i++ ) {
        if (d[i] < INF) {
            cout << d[i] << endl;
        } else {
            cout << "INF" << endl;
        }
    }

    return 0;
}

#endif