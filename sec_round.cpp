#include <bits/stdc++.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
using namespace std;

//#define TEST

#define Inline __inline__ __attribute__((always_inline))

struct X {
    char t, o;
};

#define P(T) T, '0',  T, '1', T, '2', T, '3', T, '4', T, '5', T, '6', T, '7', T, '8', T, '9'

static const X s_pairs[] = { P('0'), P('1'), P('2'), P('3'), P('4'), P('5'), P('6'), P('7'), P('8'), P('9') };

#define W(N, I) *(X*)&b[N] = s_pairs[I]

#define A(N) t = (uint64_t(1) << (32 + N / 5 * N * 53 / 16)) / uint32_t(1e##N) + 1 + N/6 - N/8, t *= u, t >>= N / 5 * N * 53 / 16, t += N / 6 * 4, W(0, t >> 32)

#define S(N) b[N] = char(uint64_t(10) * uint32_t(t) >> 32) + '0'

#define D(N) t = uint64_t(100) * uint32_t(t), W(N, t >> 32)

#define L0 b[0] = char(u) + '0'

#define L1 W(0, u)

#define L2 A(1), S(2)

#define L3 A(2), D(2)

#define L4 A(3), D(2), S(4)

#define L5 A(4), D(2), D(4)

#define L6 A(5), D(2), D(4), S(6)

#define L7 A(6), D(2), D(4), D(6)

#define L8 A(7), D(2), D(4), D(6), S(8)

#define L9 A(8), D(2), D(4), D(6), D(8)

#define LN(N) (L##N, b += N + 1)

#define LZ(N) (L##N, length = N)

#define LG(F) (u<100 ? u<10 ? F(0) : F(1) : u<1000000 ? u<10000 ? u<1000 ? F(2) : F(3) : u<100000 ? F(4) : F(5) : u<100000000 ? u<10000000 ? F(6) : F(7) : u<1000000000 ? F(8) : F(9))

Inline uint8_t myItoa(uint32_t u, char* b){
    uint8_t length;
    uint64_t t;
    LG(LZ);
    return length + 1;
}

Inline uint32_t myAtoi(char* s){
    char* t = s;
    uint32_t ans = 0;
    do{
        ans = (ans<<3) + (ans<<1);
        ans += (*t - 48);
        t++;
    }while(47 < *t && *t < 58);
    return ans;
}

Inline uint64_t myAtoi1(char* s){
    char* t = s;
    uint64_t ans = 0;
    do{
        ans = (ans<<3) + (ans<<1);
        ans += (*t - 48);
        t++;
    }while(47 < *t && *t < 58);
    ans *= 100;
    if(*t == '.'){
        t++;
        ans += (*t-48) * 10;
        t++;
        if(47 < *t && *t < 58){
            ans += (*t-48);
        }
    }
    return ans;
}

const uint8_t THREAD_NUM = 8;
uint32_t nodeCnt;
char* buf;
vector<uint32_t> unhashID;
vector<uint32_t> nodeData[THREAD_NUM];
vector<uint64_t> edgeData[THREAD_NUM];
struct Bi{
    uint32_t tid;
    uint32_t l, r;
};
vector<vector<tuple<uint32_t, uint64_t>>> Graph;
vector<vector<tuple<uint32_t, uint64_t>>> VerGraph;


Inline void readTestData(Bi& bi){
    char* bg = buf + bi.l;
    const char* ed = buf + bi.r;
    char tmp[15];
    nodeData[bi.tid].reserve(1000000);
    edgeData[bi.tid].reserve(500000);
    do{
        for(uint8_t i = 0; i < 3; ++i){
            uint8_t index = 0;
            while(*(bg+index) != ',' && *(bg+index) != '\n') index++;
            memset(tmp, ' ', 15);
            memcpy(tmp, bg, index);
            i < 2 ? nodeData[bi.tid].emplace_back(myAtoi(tmp)) : edgeData[bi.tid].emplace_back(myAtoi1(tmp));
            bg += (index + 1);
        }
    }while(bg < ed);
}

Inline void read(const string& fileName){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    int fd = open(fileName.c_str(), O_RDONLY);
    uint32_t length = lseek(fd, 0, SEEK_END);
    buf = (char*)mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
    thread threads[THREAD_NUM];
    Bi bi[THREAD_NUM];
    uint32_t pre_p = 0;
    for(uint8_t i = 0; i < THREAD_NUM; ++i){
        uint32_t p = (i + 1) * length / THREAD_NUM;
        while (buf[p - 1] != '\n') --p;
        bi[i].l = pre_p;
        bi[i].r = p;
        pre_p = p;
        bi[i].tid = i;
        threads[i] = thread(readTestData, ref(bi[i]));
    }
    for(uint8_t i = 0; i < THREAD_NUM; ++i){
        threads[i].join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "read test data time: " << dr_ms  << " ms\n";
#endif
}

Inline void buildGraph(){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    vector<uint32_t> nodeD;
    vector<uint64_t> edgeD;
    for(uint8_t i = 0; i < THREAD_NUM; ++i){
        nodeD.insert(nodeD.end(), nodeData[i].begin(), nodeData[i].end());
        edgeD.insert(edgeD.end(), edgeData[i].begin(), edgeData[i].end());
    }
    vector<uint32_t> tmp = nodeD;
    sort(tmp.begin(), tmp.end());
    tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
    nodeCnt = tmp.size();
    Graph = vector<vector<tuple<uint32_t,uint64_t>>>(nodeCnt);
    VerGraph = vector<vector<tuple<uint32_t,uint64_t>>>(nodeCnt);
    unhashID = tmp;
    unordered_map<uint32_t, uint32_t> hashID;
    for(uint32_t i = 0; i < nodeCnt; ++i)
        hashID[tmp[i]] = i;
    uint32_t send, send1, recv, recv1;
    for(uint32_t i = 0, j = 0; i < nodeD.size(); i+=4, j+=2){
        send = hashID[nodeD[i]];
        recv = hashID[nodeD[i+1]];
        Graph[send].emplace_back(make_tuple(recv, edgeD[j]));
        VerGraph[recv].emplace_back(make_tuple(send, edgeD[j]));
        send1 = hashID[nodeD[i+2]];
        recv1 = hashID[nodeD[i+3]];
        Graph[send1].emplace_back(make_tuple(recv1, edgeD[j+1]));
        VerGraph[recv1].emplace_back(make_tuple(send1, edgeD[j+1]));
    }
    for(uint32_t i = 0; i < nodeCnt; ++i){
        if(Graph[i].size() > 1) sort(Graph[i].begin(), Graph[i].end());
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "build graph time: " << dr_ms  << " ms\n";
#endif
}

atomic<uint32_t> thNode(0);
const uint32_t MAX_ID = 4000000;
uint8_t tidBelong[MAX_ID];
bool checkCycle[6][MAX_ID]{false};
uint32_t beginPos[6][MAX_ID];
uint32_t len[6][MAX_ID];

struct RES{
    uint32_t resNum = 0;
    uint32_t index[6];
    char* Res0 = new char[33 * 1000000];
    char* Res1 = new char[44 * 1000000];
    char* Res2 = new char[55 * 1000000];
    char* Res3 = new char[66 * 2000000];
    char* Res4 = new char[77 * 3000000];
    char* Res5 = new char[88 * 5000000];
    char* Res[6] = {Res0, Res1, Res2, Res3, Res4, Res5};
}res[THREAD_NUM];

Inline bool judgeW(const uint64_t a, const uint64_t b){
    return a <= ((b<<2) + b) && b <= ((a<<1) + a);
}

Inline void convert(const uint8_t tid, const uint8_t cycLen, uint32_t* nodeVec){
    char* p = res[tid].Res[cycLen] + res[tid].index[cycLen];
    for(uint8_t k = 0; k < cycLen+2; ++k){
        uint8_t l = myItoa(unhashID[nodeVec[k]], p);
        p += l;
        res[tid].index[cycLen] += l;
        *(p++) = ',';
        res[tid].index[cycLen]++;
    }
    uint8_t l = myItoa(unhashID[nodeVec[cycLen+2]], p);
    p += l;
    res[tid].index[cycLen] += l;
    *(p) = '\n';
    res[tid].index[cycLen]++;
}

void dfs(const uint32_t node, const uint8_t tid, uint32_t* nodeVec, bool* visited, uint8_t* distc){
    nodeVec[0] = node;
    visited[node] = true;
    tidBelong[node] = tid;
    for(uint16_t i = 0; i < Graph[node].size(); ++i){
        uint32_t ii = get<0>(Graph[node][i]);
        if(ii <= node) continue;
        else{
            nodeVec[1] = ii;
            visited[ii] = true;
            for(uint16_t j = 0; j < Graph[ii].size(); ++j){
                uint32_t jj = get<0>(Graph[ii][j]);
                if(jj <= node) continue;
                if(!judgeW(get<1>(Graph[node][i]), get<1>(Graph[ii][j]))) continue;
                if(!visited[jj]){
                    nodeVec[2] = jj;
                    visited[jj] = true;
                    for(uint16_t k = 0; k < Graph[jj].size(); ++k){
                        uint32_t kk = get<0>(Graph[jj][k]);
                        if(kk < node) continue;
                        if(!judgeW(get<1>(Graph[ii][j]), get<1>(Graph[jj][k]))) continue;
                        if(kk == node){
                            if(judgeW(get<1>(Graph[jj][k]), get<1>(Graph[node][i]))){
                                if(!checkCycle[0][node]) {
                                    checkCycle[0][node] = true;
                                    beginPos[0][node] = res[tid].index[0];
                                }
                                convert(tid, 0, nodeVec);
                                res[tid].resNum++;
                            }
                            continue;
                        }
                        if(!visited[kk]){
                            nodeVec[3] = kk;
                            visited[kk] = true;
                            for(uint16_t l = 0; l < Graph[kk].size(); ++l){
                                uint32_t ll = get<0>(Graph[kk][l]);
                                if(ll < node) continue;
                                if(distc[ll] == 0) continue;
                                if(!judgeW(get<1>(Graph[jj][k]), get<1>(Graph[kk][l]))) continue;
                                if(ll == node){
                                    if(judgeW(get<1>(Graph[kk][l]), get<1>(Graph[node][i]))){
                                        if(!checkCycle[1][node]) {
                                            checkCycle[1][node] = true;
                                            beginPos[1][node] = res[tid].index[1];
                                        }
                                        convert(tid, 1, nodeVec);
                                        res[tid].resNum++;
                                    }
                                    continue;
                                }
                                if(!visited[ll]){
                                    nodeVec[4] = ll;
                                    visited[ll] = true;
                                    for(uint16_t z = 0; z < Graph[ll].size(); ++z){
                                        uint32_t zz = get<0>(Graph[ll][z]);
                                        if(zz < node) continue;
                                        if(distc[zz] == 0 || distc[zz] == 4) continue;
                                        if(!judgeW(get<1>(Graph[kk][l]), get<1>(Graph[ll][z]))) continue;
                                        if(zz == node){
                                            if(judgeW(get<1>(Graph[ll][z]), get<1>(Graph[node][i]))){
                                                if(!checkCycle[2][node]) {
                                                    checkCycle[2][node] = true;
                                                    beginPos[2][node] = res[tid].index[2];
                                                }
                                                convert(tid, 2, nodeVec);
                                                res[tid].resNum++;
                                            }
                                            continue;
                                        }
                                        if(!visited[zz]){
                                            nodeVec[5] = zz;
                                            visited[zz] = true;
                                            for(uint16_t x = 0; x < Graph[zz].size(); ++x){
                                                uint32_t xx = get<0>(Graph[zz][x]);
                                                if(xx < node) continue;
                                                if(distc[xx] == 0 || distc[xx] == 3 || distc[xx] == 4) continue;
                                                if(!judgeW(get<1>(Graph[ll][z]), get<1>(Graph[zz][x]))) continue;
                                                if(xx == node) {
                                                    if (judgeW(get<1>(Graph[zz][x]), get<1>(Graph[node][i]))) {
                                                        if(!checkCycle[3][node]) {
                                                            checkCycle[3][node] = true;
                                                            beginPos[3][node] = res[tid].index[3];
                                                        }
                                                        convert(tid, 3, nodeVec);
                                                        res[tid].resNum++;
                                                    }
                                                    continue;
                                                }
                                                if(!visited[xx]){
                                                    nodeVec[6] = xx;
                                                    visited[xx] = true;
                                                    for(uint16_t c = 0; c < Graph[xx].size(); ++c){
                                                        uint32_t cc = get<0>(Graph[xx][c]);
                                                        if(cc < node) continue;
                                                        if(distc[cc] == 0 || distc[cc] == 2 || distc[cc] == 3 || distc[cc] == 4) continue;
                                                        if(!judgeW(get<1>(Graph[zz][x]), get<1>(Graph[xx][c]))) continue;
                                                        if(cc == node) {
                                                            if (judgeW(get<1>(Graph[xx][c]), get<1>(Graph[node][i]))) {
                                                                if(!checkCycle[4][node]) {
                                                                    checkCycle[4][node] = true;
                                                                    beginPos[4][node] = res[tid].index[4];
                                                                }
                                                                convert(tid, 4, nodeVec);
                                                                res[tid].resNum++;
                                                            }
                                                            continue;
                                                        }
                                                        if(!visited[cc]){
                                                            nodeVec[7] = cc;
                                                            for(uint16_t v = 0; ; ++v){
                                                                if(get<0>(Graph[cc][v]) < node) continue;
                                                                if(!judgeW(get<1>(Graph[xx][c]), get<1>(Graph[cc][v]))) break;
                                                                if(judgeW(get<1>(Graph[cc][v]), get<1>(Graph[node][i]))){
                                                                    if(!checkCycle[5][node]) {
                                                                        checkCycle[5][node] = true;
                                                                        beginPos[5][node] = res[tid].index[5];
                                                                    }
                                                                    convert(tid, 5, nodeVec);
                                                                    res[tid].resNum++;
                                                                }
                                                                break;
                                                            }
                                                        }
                                                    }
                                                    visited[xx] = false;
                                                }
                                            }
                                            visited[zz] = false;
                                        }
                                    }
                                    visited[ll] = false;
                                }
                            }
                            visited[kk] = false;
                        }
                    }
                    visited[jj] = false;
                }
            }
            visited[ii] = false;
        }
    }
    visited[node] = false;
    for(uint8_t i = 0; i < 6; ++i){
        len[i][node] = res[tid].index[i] - beginPos[i][node];
    }
}

Inline void getThreadResult(const uint8_t tid){
    uint32_t nodeVec[8];
    bool visited[nodeCnt]{false};
    uint8_t distc[nodeCnt];
    uint32_t* recov = new uint32_t[nodeCnt];
    uint32_t indexOfRecov = 0;
    while(1){
        uint32_t tmp = thNode++;
        if(tmp >= nodeCnt) break;
        if(Graph[tmp].size() != 0){
            distc[tmp] = 5;
            recov[indexOfRecov++] = tmp;
            for(uint16_t j = 0; j < VerGraph[tmp].size(); ++j){
                uint32_t jj = get<0>(VerGraph[tmp][j]);
                if(jj <= tmp) continue;
                else{
                    if(distc[jj] == 0) {
                        distc[jj] = 1;
                        recov[indexOfRecov++] = jj;
                    }
                    else if(distc[jj] >= 2) distc[jj] = 1;
                    for(uint16_t k = 0; k < VerGraph[jj].size(); ++k){
                        uint32_t kk = get<0>(VerGraph[jj][k]);
                        if(kk <= tmp) continue;
                        if(!judgeW(get<1>(VerGraph[jj][k]), get<1>(VerGraph[tmp][j]))) continue;
                        else{
                            if(distc[kk] == 0){
                                distc[kk] = 2;
                                recov[indexOfRecov++] = kk;
                            }
                            else if(distc[kk] >= 3) distc[kk] = 2;
                            for(uint16_t x = 0; x < VerGraph[kk].size(); ++x){
                                uint32_t xx = get<0>(VerGraph[kk][x]);
                                if(xx <= tmp || xx == jj) continue;
                                if(!judgeW(get<1>(VerGraph[kk][x]), get<1>(VerGraph[jj][k]))) continue;
                                if(distc[xx] == 0) {
                                    distc[xx] = 3;
                                    recov[indexOfRecov++] = xx;
                                }
                                else if(distc[xx] == 4) distc[xx] = 3;
                                for(uint16_t c = 0; c < VerGraph[xx].size(); ++c){
                                    uint32_t cc = get<0>(VerGraph[xx][c]);
                                    if(cc <= tmp || cc == kk || cc == jj) continue;
                                    if(!judgeW(get<1>(VerGraph[xx][c]), get<1>(VerGraph[kk][x]))) continue;
                                    if(distc[cc] == 0){
                                        distc[cc] = 4;
                                        recov[indexOfRecov++] = cc;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            dfs(tmp, tid, nodeVec, visited, distc);
            for(uint32_t j = 0; j < indexOfRecov; ++j) distc[recov[j]] = 0;
            indexOfRecov = 0;
        }
    }
}

Inline void getResult(){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    thread threads[THREAD_NUM];
    for(uint8_t i = 0; i < THREAD_NUM; ++i){
        threads[i] = thread(getThreadResult, i);
    }
    for(uint8_t i = 0; i < THREAD_NUM; ++i){
        threads[i].join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "get result time: " << dr_ms  << " ms\n";
#endif
}


Inline void saveResult(const string& fileName){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    uint32_t totalResultNum = 0;
    for(uint32_t i = 0; i < THREAD_NUM; ++i){
        totalResultNum += res[i].resNum;
    }
#ifdef TEST
    cout << "result size: " << totalResultNum << '\n';
#endif
    FILE* fp;
    fp = fopen(fileName.c_str(), "w");
    char* title = new char[20];
    sprintf(title, "%d\n", totalResultNum);
    fwrite(title, sizeof(char), strlen(title), fp);
    for(uint8_t i = 0; i < 6; ++i){
        for(uint32_t j = 0; j < nodeCnt; ++j){
            if(!checkCycle[i][j]) continue;
            fwrite(res[tidBelong[j]].Res[i]+beginPos[i][j], sizeof(char), len[i][j], fp);
        }
    }
    fclose(fp);
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "output result time: " << dr_ms  << " ms\n";
#endif
}


int main() {
    string testFile = "/data/test_data.txt";
    string resultFile = "/projects/student/result.txt";
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    read(testFile);
    buildGraph();
    getResult();
    saveResult(resultFile);
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double, std::milli>(end - start).count();
    cout << "total cost time: " << dr_ms << " ms\n";
#endif
    return 0;
}
