
//日付
//目標時間：
//かかった時間：

//*ポイント

/*
*KPLTパート
K：
P：1列目にA,2列目にBの入力の仕方がわからず。
L：
T：


*/

//ーーー　コードスタート　ーーー

#include <iostream>
#include <string>
using namespace std;

int main(){

int N,K,P;

cin >> N >> K >> P;

int C[N],A[N][K];

for(int i=0;i<N;i++){

    for(int j=0;j<K;j++){

        cin >> C[i];
        cin >> A[i][j];
    
    }

}

int judge =0;
//全ての施策を実行後の全てのパラーメータがP以上いくなら
for (int i = 0; i < K; i++)
{
    /* code */
    int sum_K =0;
    for (int j = 0; j < N; i++)
    {
        /* code */
        sum_K += A[j][i];
    }
    
    if (sum_K >= P)//P以上か？
    {
        judge += 1;
        /* code */
    }
    
}

//全て目標達成できているなら
if(judge ==3){

}

//達成できない
else{
    cout << -1 << endl;
}


}

