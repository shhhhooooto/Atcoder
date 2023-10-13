//日付　9/25
//目標時間：15min
//かかった時間：20分かけてもできず　9/25。次回続きから。
//30minかけてかなりいい感じになったがいまだにWA 9/26

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
#include <iostream>

using namespace std;

int main(){

string S,T;

//パターン：erase eraser dream dreamer

cin >> S;

//探索効率を考える

//単語数分繰り返し
std::string::size_type i=0;//単語の区切りに使用
bool yes_flag = 0;
while(S[0+i] == 'e' || S[0+i] == 'd'){//最初の文字がeかd

    if(S.compare(0+i,5,"erase") == 0){//e→rase
        
        //cout << "OK" << endl;

        if(S.size() == i+5){// earseの次の文字が空→fin
            i += 5;
            yes_flag = 1;
            //cout << "OK1" << endl;
            break;
        }
        
        else if(S[0+i+5] == 'r'){//r
            i += 6;
            if(S.size() == i ){
                yes_flag = 1;
                //cout << "OK2" << endl;
                break;
            }
        }

        else{//また最初から
            //cout << "OK3" << endl;
            i += 5;
        }
    }

    if(S.compare(0+i,5,"dream") == 0){//d→ream
        //cout << "Good" << endl;
        if(S.size() == i+5){// dreamの次の文字が空→fin
            i += 5;
            yes_flag = 1;
            //cout << "Good1" << endl;
            break;
        }
        
        else if(S.compare(0+i,2,"er")){//r
            i += 7;
            //cout << "Good2" << endl;
            if(S.size() == i ){
                yes_flag = 1;
                //cout << "Good2-1" << endl;
                break;
            }
        }
        
        else{//また最初から
            i += 5;
            //cout << "Good3" << endl;
        }

    }
}

// cout << S.size() << endl;
// cout << i << endl;

if(yes_flag){
    cout << "Yes" << endl;
}

else{
    cout << "No" << endl;
}

}

