#include <iostream>
#include <iomanip>
#include <algorithm>
#include <climits> 
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;

enum Dir {
    EAST = 1, WEST, NORTH, SOUTH, UP, DOWN
};

class Dice {
public:
    Dice(const int label[6])  {
        _e = 3; _n = 5; _u = 1;
        for (int d = 0; d < 6; d++) _label[d] = label[d];
    }

    int GetIndex(Dir dir) const 
    {
        switch (dir) {
        case (EAST):
            return _e;
        case (WEST):
            return 7 - _e;
        case (NORTH):
            return _n;
        case (SOUTH):
            return 7 - _n;
        case (UP):
            return _u;
        case (DOWN):
            return 7 - _u;
        default:
            return -1; // does not reach here
        }
    }

    Dir GetDir(int id) const {
        if (id == _e) return EAST;
        if (id == 7 - _e) return WEST;
        if (id == _n) return NORTH;
        if (id == 7 - _n) return SOUTH;
        if (id == _u) return UP;
        return DOWN;
    }

    void RollTo(Dir dir)
    {
        int e_new = _e, n_new = _n, u_new = _u;
        switch (dir) {
        case (EAST):
            e_new = GetIndex(UP); u_new = GetIndex(WEST); break;
        case (WEST):
            e_new = GetIndex(DOWN); u_new = GetIndex(EAST); break;
        case (NORTH):
            n_new = GetIndex(UP); u_new = GetIndex(SOUTH); break;
        case (SOUTH):
            n_new = GetIndex(DOWN); u_new = GetIndex(NORTH); break;
        default: break;
        }

        // update
        _e = e_new; _n = n_new; _u = u_new;
    }

    int GetLabel(int dice) const {  // dice index = 1-6
        return _label[dice-1];
    }
    int GetLabel(Dir dir) const {
        return GetLabel(GetIndex(dir));
    }
    void SetLabel(const int label[6])
    {
        for (int d = 0; d < 6; d++) {
            _label[d] = label[d];
        }
    }
    vector<int> FindIndex(int label) const 
    {
        vector<int> v;
        for (int i = 0; i < 6; i++) {
            if (_label[i] == label) {
                v.push_back(i+1);
            }
        }
        return v;
    }

    void RollDirToUpper(Dir dir)
    {
        if (dir == EAST) RollTo(WEST);
        if (dir == WEST) RollTo(EAST);
        if (dir == NORTH) RollTo(SOUTH);
        if (dir == SOUTH) RollTo(NORTH);
        if (dir == DOWN) {
            RollTo(NORTH);
            RollTo(NORTH);
        }
    }

private:
    int _e, _n, _u;  // dice index for EAST(x), NORTH(y), UPPER(z) , =1, ..., 6
    int _label[6];   // label for each index, _label[d-1], d=1,..,6
};

bool IsSameDice(const Dice& dice1, Dice& dice2)  // note: Dice2 will be rotated, but intrinsically not changed
{
    int upL = dice1.GetLabel(UP);
    vector<int> idxs = dice2.FindIndex(upL);
    if (idxs.empty()) {
        return false;
    }

    for ( const auto& idx : idxs ) {
        // rotate dice2 to the same UP label posi
        dice2.RollDirToUpper(dice2.GetDir(idx));

        // check DOWN Label
        if (dice1.GetLabel(DOWN) != dice2.GetLabel(DOWN)) {
            continue;
        }

        // check the othr 4
        vector<int> vL1(4, -1), vL2(4, -1);
        Dir dirs[] = {EAST, NORTH, WEST, SOUTH};
        for (int i = 0; i < 4; i++) {
            vL1[i] = dice1.GetLabel(dirs[i]);
        }
        for (int i = 0; i < 4; i++) {
            vL2[i] = dice2.GetLabel(dirs[i]);
        }

        for (int i = 0; i < 4; i++ ) {
            if (i > 0) {
                rotate(vL2.begin(), vL2.begin()+1, vL2.end());
            }
            if (vL1 == vL2) {
                return true;
            }
        }
    }
    return false;
}

int main()
{
    int n;
    cin >> n;

    int L[6];
    for (int i = 0; i < 6; i++) {
        cin >> L[i];
    }
    Dice dice1(L);

    bool found_same = false;

    for (int i = 0; i < n-1; i++) {
        for (int i = 0; i < 6; i++) {
            cin >> L[i];
        }
        Dice dice2(L);

        if (!found_same && IsSameDice(dice1, dice2)) {
            found_same = true;
        } 
    }

    if (found_same == false) {
        cout << "Yes" << endl;
    } else {
        cout << "No" << endl;
    }

    return 0;
}

#if 0

int main()
{
    int L[6];
    for (int i = 0; i < 6; i++) {
        cin >> L[i];
    }
    Dice dice1(L);

    for (int i = 0; i < 6; i++) {
        cin >> L[i];
    }
    Dice dice2(L);

    if (IsSameDice(dice1, dice2)) {
        cout << "Yes" << endl;
    } else {
        cout << "No" << endl;
    }

    return 0;
}

int main()
{
    int L[6];
    for (int i = 0; i < 6; i++) {
        cin >> L[i];
    }

    Dice dice(L);

    int q = 0;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int L1, L2;
        cin >> L1 >> L2;

        Dir L1Dir = dice.GetDir(dice.FindIndex(L1));

        dice.RollDirToUpper(L1Dir);

        Dir L2Dir = dice.GetDir(dice.FindIndex(L2));

        int L3 = -1;
        if (L2Dir == EAST) {
            L3 = dice.GetLabel(NORTH);
        } else if (L2Dir == NORTH) {
            L3 = dice.GetLabel(WEST);
        } else if (L2Dir == WEST) {
            L3 = dice.GetLabel(SOUTH);
        } else if (L2Dir == SOUTH) {
            L3 = dice.GetLabel(EAST);
        } 

        cout << L3 << endl;
    }

    return 0;
}


int main()
{
    int L[6];
    for (int i = 0; i < 6; i++) {
        cin >> L[i];
    }
    string roll;
    cin >> roll;

    Dice dice(L);

    for (int i = 0; i < roll.size(); i++) {
        Dir dir = UP;
        if (roll[i] == 'E') dir = EAST;
        if (roll[i] == 'W') dir = WEST;
        if (roll[i] == 'N') dir = NORTH;
        if (roll[i] == 'S') dir = SOUTH;

        dice.RollTo(dir);
    }

    cout << dice.GetLabel(UP) << endl;

    return 0;
}


double calc_n_dist(int p, const vector<double>& x, const vector<double>& y )
{
    double D = 0;

    int n = (int)x.size();

    if (p == 0) {
        for (int i = 0; i < n; i++ ) {
            D = max(D, fabs(x[i] - y[i]));
        }
    }
    else {
        for (int i = 0; i < n; i++ ) {
            D += pow( fabs(x[i] - y[i]), p );
        }
        double pp = 1/(double)p;
        D = pow(D, pp);
    }

    return D;
}

int main()
{
    int n;
    cin >> n;
    vector<double> x(n, 0), y(n, 0);

    for (int i = 0; i < n; i++) {
        cin >> x[i];
    }
    for (int i = 0; i < n; i++) {
        cin >> y[i];
    }

    cout << fixed << setprecision(8);
    cout << calc_n_dist( 1, x, y ) << endl;
    cout << calc_n_dist( 2, x, y ) << endl;
    cout << calc_n_dist( 3, x, y ) << endl;
    cout << calc_n_dist( 0, x, y ) << endl;

    return 0;
}

int main()
{
    for (int dat = 0; dat < 1000; dat++ ) {
        int n;
        cin >> n;

        if (n == 0) return 0;

        double avl = 0, avl2 = 0;
        for (int i = 0; i < n; i++) {
            double s;
            cin >> s;
            avl += s;
            avl2 += s * s;
        }

        avl /= (double)n;
        avl2 /= (double)n;

        cout << fixed << setprecision(8) << sqrt(avl2 - avl * avl) << endl;

    }

    return 0;
}

int main()
{
    const double PI=3.14159265358979323846;

    int ia, ib, iC;
    cin >> ia >> ib >> iC;

    double a = (double)ia;
    double b = (double)ib;
    double C = (iC / (double)180 ) * PI;
    double sin_C = sin(C);
    double cos_C = cos(C);

    double S = 0.5 * a * b * sin_C;
    double c = sqrt(a*a + b*b - 2.*a*b*cos_C);
    double L = a + b + c;
    double h = 2 * S / a;

    cout << fixed << setprecision(8) << S << endl;
    cout << L << endl;
    cout << h << endl;

    return 0;
}


int main()
{
    double x1, y1, x2, y2;
    cin >> x1 >> y1 >> x2 >> y2;

    cout << fixed << setprecision(8) << sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)) << endl;

    return 0;
}

int main()
{
    string str;
    cin >> str;

    int q;
    cin >> q;

    for (int i = 0; i < q; i++) {
        string command;
        int a, b;
        cin >> command >> a >> b;
        if (command == "print") {
            cout << str.substr(a, b - a + 1) << endl;
        }
        if (command == "reverse") {
            string ssub = str.substr(a, b - a + 1);
            reverse(ssub.begin(), ssub.end());
            str.replace( a, b - a + 1, ssub );
        }
        if (command == "replace") {
            string p;
            cin >> p;
            str.replace( a, b - a + 1, p );
        }
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    int t = 0, h = 0;
    for (int i = 0; i < n; i++) {
        string s1, s2;
        cin >> s1 >> s2;

        if (s1 > s2) {
            t += 3;
        }
        else if (s1 < s2) {
            h += 3;
        }
        else {
            t += 1;
            h += 1;
        }
    }

    cout << t << " " << h << endl;

    return 0;
}

int main()
{
    for (int ii = 0; ii < 10; ii++) {

        char cards[210] = {0};
        cin >> cards;

        if (cards[0] == '-') {
            return 0;
        }

        int m;
        cin >> m; 

        int n = (int)strlen(cards);

        char buf[210] = {0};
        for (int i = 0; i < m; i++) {
            int h;
            cin >> h;
            for (int j = 0; j < h; j++) {
                buf[n-h+j] = cards[j];
            }
            for (int j = 0; j < n-h; j++) {
                buf[j] = cards[h+j];
            }
            swap(cards, buf);
        }

        cout << cards << endl;
    }

    return 0;
}

int s_stricmp(const char *s1, const char *s2)
{
  register const unsigned char *ss1, *ss2;
  for (ss1 = (const unsigned char*)s1, ss2 = (const unsigned char*)s2;
       *ss1 != '\0' && tolower(*ss1) == tolower(*ss2) ;
       ss1++, ss2++)
    ;
  return *ss1 - *ss2;
}

int main()
{
    string target;
    cin >> target;

    int count = 0;
    while (1) {
        string w;
        cin >> w;

        if (w == "END_OF_TEXT")
            break;

        if (s_stricmp(target.c_str(), w.c_str()) == 0) {
            count++;
        }
    }

    cout << count << endl;

    return 0;
}


int main()
{
    string s, p;
    getline(cin, s);
    getline(cin, p);

    size_t ns = s.size();
    size_t np = p.size();

    for (size_t i = 0; i < ns; i++) {
        rotate(s.begin(), s.begin()+1, s.end());
        if ( s.substr(0, np) == p ) {
            cout << "Yes" << endl;
            return 0;
        }
    }

    cout << "No" << endl;

    return 0;
}

int main()
{
    map<char, int>  charcount;
    for ( char cc = 'a'; cc <= 'z'; cc++ ) {
        charcount[cc] = 0;
    }

    char cstr[1201];
    while ( cin >> cstr ) {
        size_t n = strlen(cstr);
        for ( size_t i = 0; i < n; i++ ) {
            char c = cstr[i];
            if (isalpha(c)) {
                c = tolower(c);
                charcount[c] ++;
            }
        }
    }

    for (auto& i : charcount) {
        cout << i.first << " : " << i.second << endl;
    }

    return 0;
}

int main()
{
    while (1) {
        string  str;
        getline(cin, str);

        int x = atoi(str.c_str());
        if ( x == 0 ) {
            return 0;
        }

        int sum = 0;
        size_t n = str.size();
        for ( size_t i = 0; i < n; i++ ) {
            char c = str[i];
            if (isdigit(c)) {
                sum += ( c - '0');
            }
        }
        cout << sum << endl;
    }

    return 0;
}

int main()
{
    string  str;
    getline(cin, str);

    const char* cstr = str.c_str();

    string str2;
    size_t n = str.size();
    for ( size_t i = 0; i < n; i++ ) {
        char c = cstr[i];
        if (isalpha(c)) {
            if (islower(c)) {
                c = toupper(c);
            }
            else if (isupper(c)) {
                c = tolower(c);
            }
        }
        str2.push_back(c);
    }
    cout << str2 << endl;

    return 0;
}

int main()
{
    int n, m, l;
    cin >> n >> m >> l;

    vector<vector<long long>> A(n, vector<long long>(m, 0));
    vector<vector<long long>> B(m, vector<long long>(l, 0));
    vector<vector<long long>> C(n, vector<long long>(l, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> A[i][j];
        }
    }
    for (int j = 0; j < m; j++) {
        for (int k = 0; k < l; k++) {
            cin >> B[j][k];
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            for (int k = 0; k < m; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            cout << C[i][j];
            if ( j < l-1 ) cout << " ";
        }
        cout << endl;
    }

    return 0;
}


int main()
{
    int n, m;
    cin >> n >> m;

    vector<vector<int>> A(n+1, vector<int>(m+1, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> A[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A[i][m] += A[i][j];
        }
    }
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            A[n][j] += A[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        A[n][m] += A[i][m];
    }

    for (int i = 0; i < n+1; i++) {
        for (int j = 0; j < m+1; j++) {
            cout << A[i][j];
            if ( j < m) cout << " ";
        }
        cout << endl;
    }

    return 0;
}

void search_sum(int m, int max_selected, int n, int target, int& count)
{
    if (m == 1) {
        if (max_selected < target && target <= n) {
            count++;
        }
        return;
    }
    for (int i = max_selected+1; i <= n; i++) {
        search_sum(m-1, i, n, target - i, count);
    }
}

int main()
{
    while (1) {
        int n, x;
        cin >> n >> x;
        if (n == 0 && x == 0 )
            return 0;

        int count = 0;
        search_sum(3, 0, n, x, count);

        cout << count << endl;
    }

    return 0;
}


int main()
{
    while (1) {
        int m, f, r;
        cin >> m >> f >> r;

        if (m < 0 && f < 0 && r < 0)
            return 0;

        char grade = '?';

        if (m < 0 || f < 0) {
            grade = 'F';
        }
        else if (m + f >= 80) {
            grade ='A';
        }
        else if (m + f >= 65) {
            grade = 'B';
        }
        else if (m + f >= 50) {
            grade = 'C';
        }
        else if (m + f >= 30 ) {
            grade = ( r >= 50 )? 'C' : 'D';
        }
        else {
            grade = 'F';
        }

        cout << grade << endl;
    }

    return 0;
}


int main()
{
    int n, m;
    cin >> n >> m;

    vector<int> A(m*n, 0);
    vector<int> b(m, 0), c(n, 0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> A[m*i + j];
        }
    }
    for (int j = 0; j < m; j++) {
        cin >> b[j];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            c[i] += A[m*i + j] * b[j];
        }
    }

    for (int i = 0; i < n; i++) {
        cout << c[i] << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    int bfr[4][3][10] = {0};

    for (int i = 0; i < n; i++) {
        int b, f, r, v;
        cin >> b >> f >> r >> v;
        bfr[b-1][f-1][r-1] += v;
    }
    for (int b = 0; b < 4; b++) {
        for (int f = 0; f < 3; f++) {
            for (int r = 0; r < 10; r++) {
                cout << " ";
                cout << bfr[b][f][r];
            }
            cout << endl;
        }
        if (b < 3) cout << "####################" << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    vector<bool> c(52, false);  // S1-13, H1-13, C1-13, D1-13
    for (int i = 0; i < n; i++) {
        char m;
        int  n;
        cin >> m >> n;
        int card = n-1;
        switch (m) {
        case 'D':
            card += 13;
        case 'C':
            card += 13;
        case 'H':
            card += 13;
        }
        c[card] = true;
    }
    for (int i = 0; i < 52; i++) {
        if (c[i]) continue;
        int im = i / 13;
        char m;
        if ( im == 0 ) m = 'S';
        if ( im == 1 ) m = 'H';
        if ( im == 2 ) m = 'C';
        if ( im == 3 ) m = 'D';

        cout << m << " " << i - 13*im + 1 << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    vector<int> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }
    for (int i = n-1; i >= 0; i--) {
        cout << v[i];
        if (i > 0) cout << " ";
    }
    cout << endl;

    return 0;
}

void call(int n){
    for (int i = 1; i <= n; i++ ) {
        int x = i;
        if ( x % 3 == 0 ){
            cout << " " << i;
            continue;
        }
        while ( x > 0) {
            if ( x % 10 == 3 ){
                cout << " " << i;
                break;
            }
            x /= 10;
        }
    }
    cout << endl;
}

int main()
{
    int n;
    cin >> n;
    call(n);

    return 0;
}

int main()
{
    while (1) {
        int h, w;
        cin >> h >> w;
        if (h<=0 && w<=0) break;

        for ( int i = 0; i < h; i++) {
            for ( int j = 0; j < w; j++) {
                if ((i + j)%2 == 0)
                    cout << "#";
                else
                    cout << ".";
            }
            cout << endl;
        }
        cout << endl;
    }

    return 0;
}

int main()
{
    while (1) {
        int h, w;
        cin >> h >> w;
        if (h<=0 && w<=0) break;

        for ( int i = 0; i < h; i++) {
            for ( int j = 0; j < w; j++) {
                if (i == 0 || i == h-1 || j == 0 || j == w-1)
                    cout << "#";
                else
                    cout << ".";
            }
            cout << endl;
        }
        cout << endl;
    }

    return 0;
}

int main()
{
    while (1) {
        int h, w;
        cin >> h >> w;
        if (h<=0 && w<=0) break;

        for ( int i = 0; i < h; i++) {
            for ( int j = 0; j < w; j++) {
                cout << "#";
            }
            cout << endl;
        }
        cout << endl;
    }

    return 0;
}

int main()
{
    int n;
    long long m = LLONG_MAX, M = LLONG_MIN;
    long long s = 0;
    cin >> n;
    for ( int i = 0; i < n; i++) {
        long long a;
        cin >> a;
        m = min(m, a);
        M = max(M, a);
        s += a;
    }

    cout << m << " " << M << " " << s << endl;

    return 0;
}

int main()
{
    int a, b;
    char op;
    while (1) {
        cin >> a >> op >> b;

        int res = 0;
        switch(op) {
        case '+':
            res = a + b; break;
        case '-':
            res = a - b; break;
        case '*':
            res = a*b; break;
        case '/':
            res = a/b; break;
        default:
            return 0;
        }

        cout << res << endl;
    }

    return 0;
}

int main()
{
    const double PI=3.14159265358979323846;

    double r;
    cin >> r;

    cout << fixed << setprecision(8) << PI*r*r << " " << 2*PI*r << endl;

    return 0;
}

int main()
{
    int a, b;
    cin >> a >> b;

    cout << a/b << " " << a - b*(a/b) << " " << fixed << setprecision(10) << (double)a/(double)b <<endl;

    return 0;
}

int main()
{
    int count = 0;
    int a, b, c;
    cin >> a >> b >> c;
    for ( int i = a; i <= b; i++) {
        if ( c % i == 0 ) {
            count++;
        }
    }

    cout << count << endl;

    return 0;
}

int main()
{
    while (1) {
        int x, y;
        cin >> x >> y;
        if (x == 0 && y == 0) break;
        cout << min(x, y) << " " << max(x, y) << endl;
    }

    return 0;
}


int main()
{
    int x[10001] = {0};
    for ( int i = 0; i < 10001; i++ ) {
        cin >> x[i];
        if (x[i] == 0 ) break;
    }

    for ( int i = 0; i < 10001; i++ ) {
        if (x[i] == 0) break;
        cout << "Case " << i+1 << ": " << x[i] << endl;
    }

    return 0;
}

int main()
{
    for ( int i = 0; i < 1000; i++ )
        cout << "Hello World" << endl;

    return 0;
}

int main()
{
    int w, h, x, y, r;
    cin >> w >> h >> x >> y >> r;

    bool yes = true;
    if (w < 2*r) {
        yes = false;
    }
    else if (h < 2*r) {
        yes = false;
    }
    else {
        if (x < r || x > w-r) {
            yes = false;
        }
        if (y < r || y > h-r) {
            yes = false;
        }
    }

    if (yes) {
        cout << "Yes" << endl;
    } else {
        cout << "No" << endl;
    }

    return 0;
}

int main()
{
    int a[3];
    for ( int i = 0; i < 3; i++)
        cin >> a[i];

    sort(a+0, a+3);

    cout << a[0] << " " << a[1] << " " << a[2] << endl;

    return 0;
}

int main()
{
    int a, b, c;
    cin >> a >> b >> c;

    if ( a < b && b < c )
        cout << "Yes" << endl;
    else
        cout << "No" << endl;

    return 0;
}


int main()
{
    int a, b;
    cin >> a >> b;

    if ( a < b )
        cout << "a < b" << endl;
    else if ( a > b )
        cout << "a > b" << endl;
    else
        cout << "a == b" << endl;

    return 0;
}


int main()
{
    int h, m, s;
    cin >> s;
    h = s / 3600;
    s = s - 3600 * h;
    m = s / 60;
    s = s - 60 * m;
    
    cout << h << ":" << m << ":" << s << endl;

    return 0;
}

int main()
{
    int a, b;
    cin >> a >> b;
    cout << a * b << " " << 2 * (a + b) << endl;

    return 0;
}

int main()
{
    int x;
    cin >> x;
    cout << x * x * x << endl;

    return 0;
}

int main()
{
    cout << "Hello World" << endl;

    return 0;
}

#endif

