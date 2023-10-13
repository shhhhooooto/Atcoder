//日付
//目標時間：30min
//かかった時間：40min

//*ポイント：計算量を考慮したプログラミング

/*
*KPLTパート
K：breakを使って動作時間を早くした。Atcoderでも動作時間over？はあるみたい
P：forループの上限値をNではなくnum_man,num_gosen,num_senとしてしまった。
そもそもこの変数やman,gosen,senの変数いらなかった。
!計算量を考慮する必要がある

L：1秒間にできるforの計算は10^8
お金の計算は0〜10枚みたいな時は、forを0<10+1にする必要ある
T：余計な変数を用意しない。ただ用意せずにややこしくなるのは避ける。いらなかったらあとで修正するスタンス
*/

#include <iostream>
using namespace std;

int main(){

int man,gosen,sen;
int N,y;
int num_man,num_gosen,num_sen;
int man_out,gosen_out,sen_out;

man = 10000;
gosen = 5000;
sen = 1000;


cin >> N;
cin >> y;

//どれを出力しても正解なので、配列の0要素を取り出すとかでいい
bool nothing_flag = true;//一つでも合ったら0を代入

for(int i=0;i<N+1;i++){

    for(int j=0;j<N+1;j++){
    
            int k = N-i-j;

            if(y == i * man + j * gosen + k * sen&& k>=0){

                man_out = i;
                gosen_out = j;
                sen_out = k;
                
                nothing_flag = 0;

                break;
            }

    }

    if (nothing_flag==0)
    {
        break;
    }
}

if (nothing_flag==0)
{
    cout << man_out << " " << gosen_out << " " << sen_out << endl;
}

else
{
    cout << -1 << " " << -1 << " " << -1 << endl;
}



    
}