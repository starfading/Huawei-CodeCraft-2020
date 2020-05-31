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

#define LG(F) (u<100 ? u<10 ? F(0) : F(1) : u<10000 ? u<1000 ? F(2) : F(3) : F(4))


Inline int myItoa(uint32_t u, char* b){
    int length;
    uint64_t t;
    LG(LZ);
    return length + 1;
}

Inline int myAtoi(char* s){
    char* tmp = s;
    register int ans = 0;
    do{
        ans = (ans<<3) + (ans<<1);
        ans += (*tmp - 48);
        tmp++;
    }while(47 < *tmp && *tmp < 58);
    return ans;
}

int rawData[560000];
int indexOfRawData;

Inline void readTestData(const string& fileName){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    int fd = open(fileName.c_str(), O_RDONLY);
    int length = lseek(fd, 0, SEEK_END);
    char* buf = (char*)mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
    char* bg = buf;
    const char* ed = buf + length;
    register int index;
    char tmp[10];
    register int i;
    do{
        for(i = 0; i < 3; ++i){
            index = 0;
            while(*(bg+index) != ',' && *(bg+index) != '\n') index++;
            if(i < 2){
                memset(tmp, ' ', 10);
                memcpy(tmp, bg, index);
                rawData[indexOfRawData++] = myAtoi(tmp);
            }
            bg += (index + 1);
        }
    }while(bg < ed);
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "read test data time: " << dr_ms  << " ms\n";
#endif
}

const int MAX_ID_NUM = 50000;
int Graph[MAX_ID_NUM][50];
int indexOfGraph[MAX_ID_NUM];
int VerGraph[MAX_ID_NUM][50];
int indexOfVerGraph[MAX_ID_NUM];
int collectionOfGraphNode[MAX_ID_NUM];
int indexOfCOGN;

Inline void buildGraph(){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    register int i;
    for(i = 0; i < indexOfRawData; i+=4){
        if(rawData[i] < MAX_ID_NUM && rawData[i+1] < MAX_ID_NUM){
            Graph[rawData[i]][indexOfGraph[rawData[i]]++] = rawData[i+1];
            VerGraph[rawData[i+1]][indexOfVerGraph[rawData[i+1]]++] = rawData[i];
        }
        if(rawData[i+2] < MAX_ID_NUM && rawData[i+3] < MAX_ID_NUM){
            Graph[rawData[i+2]][indexOfGraph[rawData[i+2]]++] = rawData[i+3];
            VerGraph[rawData[i+3]][indexOfVerGraph[rawData[i+3]]++] = rawData[i+2];
        }
    }
    for(i = 0; i < MAX_ID_NUM; ++i){
        if(indexOfGraph[i] != 0) {
            collectionOfGraphNode[indexOfCOGN++] = i;
            if(indexOfGraph[i] > 1) sort(Graph[i], Graph[i] + indexOfGraph[i]);
        }
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "build graph time: " << dr_ms  << " ms\n";
#endif
}


const int THREAD_NUM = 8;
int Result0[THREAD_NUM][3 * 500000];
int Result1[THREAD_NUM][4 * 500000];
int Result2[THREAD_NUM][5 * 1000000];
int Result3[THREAD_NUM][6 * 2000000];
int Result4[THREAD_NUM][7 * 3000000];
int* threadResult[5][THREAD_NUM] = {{Result0[0], Result0[1], Result0[2], Result0[3],Result0[4], Result0[5], Result0[6], Result0[7]},
                                    {Result1[0], Result1[1], Result1[2], Result1[3],Result1[4], Result1[5], Result1[6], Result1[7]},
                                    {Result2[0], Result2[1], Result2[2], Result2[3],Result2[4], Result2[5], Result2[6], Result2[7]},
                                    {Result3[0], Result3[1], Result3[2], Result3[3],Result3[4], Result3[5], Result3[6], Result3[7]},
                                    {Result4[0], Result4[1], Result4[2], Result4[3],Result4[4], Result4[5], Result4[6], Result4[7]}};
int indexOfTR[5][THREAD_NUM];
int threadNodeList[THREAD_NUM][MAX_ID_NUM];
int indexOfTNL[THREAD_NUM];

void searchResult(const int& curNode, int depth, int* nodeVec, bool* visited, bool* linkable, const int& tid, vector<pair<int,int>>* Bridge){
    nodeVec[depth - 1] = curNode;
    visited[curNode] = true;
    register int index, i, j, tmp;
    for(index = 0; index < indexOfGraph[curNode]; ++index) if(Graph[curNode][index] >= nodeVec[0]) break;
    for(; index < indexOfGraph[curNode]; ++index){
        tmp = Graph[curNode][index];
        if(!visited[tmp]){
            if(linkable[tmp]){
                for(i = 0; i < Bridge[tmp].size(); ++i){
                    if(visited[Bridge[tmp][i].first] || visited[Bridge[tmp][i].second]) continue;
                    for(j = 0; j < depth; ++j) threadResult[depth][tid][indexOfTR[depth][tid]++] = nodeVec[j];
                    threadResult[depth][tid][indexOfTR[depth][tid]++] = tmp;
                    threadResult[depth][tid][indexOfTR[depth][tid]++] = Bridge[tmp][i].first;
                    threadResult[depth][tid][indexOfTR[depth][tid]++] = Bridge[tmp][i].second;
                }
            }
            if(depth < 4 && indexOfGraph[tmp] != 0) searchResult(tmp, depth + 1, nodeVec, visited, linkable, tid, Bridge);
        }
    }
    visited[curNode] = false;
}

Inline bool myCmp(const int a[2], const int b[2]){
    if(a[0] != b[0]) return a[0] < b[0];
    else return a[1] < b[1];
}

Inline bool dicsort(const pair<int,int>& a, const pair<int,int>& b){
    if(a.first != b.first) return a.first < b.first;
    else return a.second < b.second;
}

Inline void getResultByThread(const int& tid){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    int nodeVec[4];
    bool visited[MAX_ID_NUM]{false};
    bool linkable[MAX_ID_NUM]{false};
    vector<pair<int,int>>* Br = new vector<pair<int,int>>[MAX_ID_NUM];
    int recover[MAX_ID_NUM];
    int indexOfRecover = 0;
    int i, j, k, x;
    int ii, jj, kk;
    int* res3[2500];
    int indexOfRes3 = 0;
    for (i = 0; i < indexOfTNL[tid]; ++i){
        ii = threadNodeList[tid][i];
        if(indexOfGraph[ii] != 0) {
            for(j = 0; j < indexOfVerGraph[ii]; ++j){
                jj = VerGraph[ii][j];
                if(jj > ii){
                    for(k = 0; k < indexOfVerGraph[jj]; ++k){
                        kk = VerGraph[jj][k];
                        if(kk > ii){
                            for(x = 0; x < indexOfVerGraph[kk]; ++x){
                                if(VerGraph[kk][x] == ii){
                                    res3[indexOfRes3] = new int[2];
                                    res3[indexOfRes3][0] = kk;
                                    res3[indexOfRes3][1] = jj;
                                    indexOfRes3++;
                                }
                                else if(VerGraph[kk][x] > ii  && VerGraph[kk][x] != jj){
                                    if(Br[VerGraph[kk][x]].empty()){
                                        linkable[VerGraph[kk][x]] = true;
                                        recover[indexOfRecover++] = VerGraph[kk][x];
                                    }
                                    Br[VerGraph[kk][x]].emplace_back(make_pair(kk, jj));
                                }
                            }
                        }
                    }
                }
            }
            if(indexOfRes3 != 0){
                if(indexOfRes3 > 1) sort(res3, res3 + indexOfRes3, myCmp);
                for(j = 0; j < indexOfRes3; ++j){
                    threadResult[0][tid][indexOfTR[0][tid]++] = ii;
                    threadResult[0][tid][indexOfTR[0][tid]++] = res3[j][0];
                    threadResult[0][tid][indexOfTR[0][tid]++] = res3[j][1];
                }
                indexOfRes3 = 0;
            }
            if(indexOfRecover != 0){
                for(j = 0; j < indexOfRecover; ++j) {
                    if(Br[recover[j]].size() > 1) sort(Br[recover[j]].begin(), Br[recover[j]].end(), dicsort);
                }
                searchResult(ii, 1, nodeVec, visited, linkable, tid, Br);
                for(j = 0; j < indexOfRecover; ++j) {
                    Br[recover[j]].clear();
                    linkable[recover[j]] = false;
                }
                indexOfRecover = 0;
            }
        }
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "thread: " << tid << " get result time: " << dr_ms  << " ms\n";
#endif
}


Inline void getResult(){
    int i, j;
    int mu = 100;
    int para[8] = {3, 5, 8, 13, 17, 25, 38, 100};
    int pre_p = 0;
    for(i = 0; i < THREAD_NUM; ++i){
        int p = para[i] * indexOfCOGN / mu;
        for(j = pre_p; j < p; ++j){
            threadNodeList[i][indexOfTNL[i]++] = collectionOfGraphNode[j];
        }
        pre_p = p;
    }
    thread threads[THREAD_NUM];
    for(i = 0; i < THREAD_NUM; ++i){
        threads[i] = thread(getResultByThread, i);
    }
    for(i = 0; i < THREAD_NUM; ++i){
        threads[i].join();
    }
}


struct Bi{
    int cycLen;
    int left, right;
    char* ans = NULL;
    int length = 0;
};

Inline void transforRes(Bi& bi, int* Result){
    bi.ans = new char[(bi.right - bi.left) * 11];
    char* p = bi.ans;
    int l;
    for(int i = bi.left; i < bi.right; ++i){
        l = myItoa(Result[i], p);
        p += l;
        bi.length += l;
        if((i + 1) % bi.cycLen == 0){
            *(p++) = '\n';
            bi.length++;
        }
        else {
            *(p++) = ',';
            bi.length++;
        }
    }
}


Inline void saveResult(const string& fileName){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
#ifdef TEST
    int totalResultNum = 0;
    for(int i = 0; i < 5; ++i){
        for(int j = 0; j < THREAD_NUM; ++j){
            totalResultNum += indexOfTR[i][j] / (i + 3);
        }
    }
    cout << "result size: " << totalResultNum << '\n';
#endif
    FILE* fp;
    fp = fopen(fileName.c_str(), "w+");
    char tit[2] = {'1', '\n'};
    fwrite(tit, sizeof(char), 2, fp);
    const int TN_Trans = 8;
    int i, k, j;
    for(i = 0; i < 5; ++i){
        for(k = 0; k < THREAD_NUM; k++){
            if(indexOfTR[i][k] < 5000){
                char* ans = new char[indexOfTR[i][k] * 11];
                char* p = ans;
                int l;
                int length(0);
                for(int ii = 0; ii < indexOfTR[i][k]; ++ii){
                    l = myItoa(threadResult[i][k][ii], p);
                    p += l;
                    length += l;
                    if((ii + 1) % (i+3) == 0){
                        *(p++) = '\n';
                        length++;
                    }
                    else {
                        *(p++) = ',';
                        length++;
                    }
                }
                fwrite(ans, sizeof(char), length, fp);
            }
            else{
                thread threads[TN_Trans];
                Bi bi[TN_Trans];
                int pre_p = 0;
                for(j = 0; j < TN_Trans; ++j){
                    int po = (indexOfTR[i][k]/(i+3)) * (j+1)/TN_Trans * (i+3);
                    bi[j].cycLen = i + 3;
                    bi[j].left = pre_p;
                    bi[j].right = po;
                    pre_p = po;
                    threads[j] = thread(transforRes, ref(bi[j]), threadResult[i][k]);
                }
                for(j = 0; j < TN_Trans; ++j) threads[j].join();
                for(j = 0; j < TN_Trans; ++j) {
                    fwrite(bi[j].ans, sizeof(char), bi[j].length, fp);
                }
            }
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
    readTestData(testFile);
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
