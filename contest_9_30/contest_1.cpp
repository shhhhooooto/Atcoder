//日付：2023/09/30
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

int N;
string S;

cin >> N;
cin >> S;

int banme;
bool flag=0;

for(int i=0;i<S.size();i++){

    if (S[i] == 'A')
    {
        /* code */
        if (S[i+1] == 'B')
        {
            /* code */
            if (S[i+2] == 'C')
            {
                /* code */
                banme = i+1;
                flag = 1;
                break;
            }
            
        }
        
    }
    
}

if (flag)
{
    /* code */
    cout << banme << endl;
}

else
{
    cout << -1 << endl;
}







//現れるならn番目


//現れないなら-1

    
}

