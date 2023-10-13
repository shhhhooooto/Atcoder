//日付：
//目標時間：
//かかった時間：

//*ポイント

/*
*KPLTパート
K：
P：
L：
T：

全て〜なら。のやり方。
〜番目と配列の要素の指定の仕方をミスった
string は数字でも"0"で判定する

*/

//ーーー　コードスタート　ーーー

#include <iostream>
#include <string>
using namespace std;

int main(){
string S;
cin >> S;

//cout << S[1] << endl;

int all_judge =0;

for(int i=1;i<9;i++){

// cout << 2*i << endl;
// cout << S[2*i] << endl;
// cout << "" << endl;
    if (S[2*i-1] == '0')
    {
        all_judge += 1;
        /* code */
    }
}

//cout << all_judge << endl;

if (all_judge == 8)
{
    /* code */
    cout << "Yes" << endl;
}

else
{
    cout << "No" << endl;
}


    
}

