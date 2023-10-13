
//日付
//目標時間：
//かかった時間：

//*ポイント

/*
*KPLTパート
K：
P：
L：
T：


*/

//ーーー　コードスタート　ーーー

#include <iostream>
#include <string>
using namespace std;

int main(){

bool pre;
bool suf;

int N,M;
string S,T;

cin >> N >> M;
cin >> S;
cin >> T;


if (S.substr(0,N) == T.substr(0,N))
{
    /* code */
    if (S.substr(0,N) == T.substr(M-N,M))//1 && 1　→0
    {
        /* code */
        cout << 0 << endl;
    }

    else//1 && 0 →1
    {
        /* code */
        cout << 1 << endl;
    }

}

else
{
    if (S.substr(0,N) == T.substr(M-N,M))//0 && 1 →2
    {
        /* code */
        cout << 2 << endl;
    }
    
    else//0 && 0 →3
    {
        /* code */
        cout << 3 << endl;
    }
    
}








    
}

