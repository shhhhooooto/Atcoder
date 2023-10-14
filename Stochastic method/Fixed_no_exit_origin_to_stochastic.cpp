/*
JetとQGPの複数散乱を記述するプログラム。を元にしたstochastic method

*データ取得用
g++ Fixed_no_exit.cpp GaussLaguerre.cxx  Function_for_col.cpp -o col.out
clang++ -std=c++11 Fixed_no_exit.cpp GaussLaguerre.cxx  Function_for_col.cpp -o col.out

todoーーーーーーーーーーーーーーーー

todo クラス使ってみる？
    →雛形を大量生産するのに向いているので、あまり実用性はないかも。
    例）人から勇者、魔女など、基本に＋＠して作っていきたいもの
    これから使うとしても、粒子衝突→2フレーバー、3フレーバ＋リコイルとかそういうパターンを作りたい時かな？

todo C
無限までのガウスらゲール使用
ゼータ関数を使用
複数粒子が箱に入ってくることも考慮する。じゃないとPYTHIAに対応できない。ただ、QGPの多重散乱を考慮すれば、これも散乱後のタイムステップから見れば、複数のQGPパートンとジェットが入射する散乱なので、同じアルゴリズムで解決できそう

todoーーーーーーーーーーーーーーーー

*/

//ーーーーーーーーーーーーープログラム開始ーーーーーーーーーーーーー

#include <fstream>
#include <math.h>
#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "GaussLaguerre.h"
#include "Function_for_col.h"

using namespace std;

//定数
double zeta_3 = 1.20205; // アペリーの数
double h_barc = 197.327; // 1 = ℏc ＝197.327MeV*fm

//*ーーーーーーー実験の設定ーーーーーーーー
int events=10;//イベント回数

//*ーーーーーーパラメータ設定ーーーーーーー
double ini_jet_P = 5000.0;//(MeV) 1000MeV=1GeV
double ini_jet_E = abs(ini_jet_P);//(MeV)*初期のjetのエネルギー(masslessなので運動量の大きさに一致)
double sigma=0.1/pow(h_barc,2.0);//散乱断面積(fm^2):運動量が主役なのでMeVの単位に合わせる
double T = 300.0; //温度 MeV
double L = 4.0/h_barc; //(fm)流体素編の一片の長さ fm=1/200(Mev)

double time_max=1.0/h_barc;//(fm)時間発展の最大値
double ts_wid=0.1/h_barc;//(fm)タイムステップ幅
int ts_max = floor(time_max/ts_wid);//タイムステップの数（tsループの回数使用）floorで切り捨て
int d_g = 16,d_q = 36; //?自由度　gluon(boson):16(8(カラー)*2(スピン？？？ スピン1の粒子なのに？))   quark (fermion):36(3(カラー)*2(スピン)*3(フレーバー)*2(クォーク/反クォーク))
int num_type_max =2;//区別する粒子の種類(0:bose、1:ferimion)

//チェック機構のo/off
bool conserve_check = 1;//保存量確認

//積分した関数F(p) ※0から∞の方に変更したい。
//積分範囲 xmin ~ xmax =F(xmax) ~ F(xmin)
double xmin = 0.0001;//どこまで積分するか
double xmax = 30*T;//大きい数字に変数変換分のTをかけておく

//ーーーーーーーーーーーーーーーーグローバル変数ーーーーーーーーーーーーーーーー
//衝突に関する変数
double d_rel_sqr=sigma/M_PI;//衝突径数の最大値（これより小さければ衝突）
int num_of_jet;//ジェットの個数

//保存量の確認
int check_all_col=0;//保存○の衝突数をカウント
int col_times_total=0;//全衝突回数をカウント→保存量の確認の際の分母に使用
int check_all_ev=0;//保存○のイベント数をカウント
double Etot_recoiled_QGP_allev = 0.0;//全イベントの差っ引かれたエネルギーの総計。最後の平均を出す時に使用

//データの取得等
int col_times_total_1ev;//1evあたりの合計衝突回数
int num_of_jet_fin;//終了時のnum_of_jetを入れて最後に表示
int num_of_out_jet_fin;//実験終了時に素片の外に出ているジェット

//Lorents transformation
vector<double> Lorentz_trans(2);//Lorentz transformation
vector<double> Lorentz_trans_2(2);//inverse transformation

// ? ここで宣言する必要ある？
vector<double> p_QGP_lab_bf(4),x_QGP_lab_bf(4);

//流体素片（箱）の中の全粒子数
//{L(fm)*T(MeV) * 1 / h_barc (=1/197.327fm*MeV)} =L*T/197.327 [無次元]
double total_mean_num_boson = pow(L * T , 3.0) * d_g * zeta_3  /pow(M_PI, 2.0) ;//boson
double total_mean_num_fermi = 3.0/4.0 * pow(L * T , 3.0) * d_q * zeta_3  /pow(M_PI, 2.0) ;//fermion
double N_total_mean = total_mean_num_boson + total_mean_num_fermi;
double E_mean = (d_g+(7.0/8.0)*d_q) * pow(L,3.0) * (pow(M_PI, 2.0)/30.0) * pow(T,4.0);//QGPのエネルギー

//*ーーーーーー分布関数(一部範囲をmain関数内で宣言)ーーーーーー
int seed=1;
std::mt19937 mt(seed);//メルセンヌツイスタを使用
//一様分布
std::uniform_real_distribution<double> p_dist_bose;
std::uniform_real_distribution<double> p_dist_fermi;
std::uniform_real_distribution<double> dis_theta(-1,1);//-1≤cosθ≤1
std::uniform_real_distribution<double> dis_phi(0,2*M_PI);//0≤φ"≤"2π
//ポアソン分布
std::uniform_real_distribution<double> distr(-L/2.0,L/2.0);//座標
std::poisson_distribution<> dist_bose(total_mean_num_boson);//boson運動量
std::poisson_distribution<> dist_fermi(total_mean_num_fermi);//fermion運動量

//*イベントごとに初期化する配列と変数
int num_of_col_ts_1ev,num_of_no_col_ts_1ev,check_1ev;
vector<int> col_ts_ID;//衝突が起きるtsを記録しておく

//*tsごとに初期化する配列
vector<vector<double> > all_QGP_x_1ts(0,vector<double> (4));//resize時に自動で初期化される。
vector<vector<double> > all_QGP_p_1ts(0,vector<double> (4));
vector<vector<double> > coll_pairs_w(0,vector<double>(4));//被りを弾く前のペアを記録する配列
vector<vector<double> > coll_pairs_1on1(0,vector<double>(4));//弾いた配列＝衝突させるペアを記録
int coll_times_w_in_ts,coll_times_in_ts,N_tot;
bool no_col_ts;
//衝突回数の次元(coll_ID),衝突時刻+衝突したペアのIDを記録する次元(tc,j,Q,ID)=(4)
vector<double> coll_jet_ID,no_coll_jet_ID;//衝突が起こる/らないjetのIDの配列、粒子数はresizeで増やすので初期サイズ1にしとく


//*ーーーーーーーーーーーーー関数ーーーーーーーーーーーーー
//割り当てる関数
void QGP_assign(vector<vector<double> > & QGP_x,vector<vector<double> > & QGP_p);
//衝突判定を行う関数
void col_judge(vector<vector<double> >& col_pair_W, vector<vector<vector<double> > > & all_jx_1ev,vector<vector<vector<double> > > & all_jp_1ev,
    int ts_lb);
//被りを除去する関数
void w_erase(vector<vector<double> >& coll_pairs_w,vector<vector<double> >& coll_pairs_1on1);
//衝突後の情報を更新する関数
void update_af_col(int i,int ts_lab,vector<double> &x_jet_lab_bf,vector<double> &p_jet_lab_bf,vector<double> &x_jet_lab_af,vector<double> &p_jet_lab_af,
vector<double> &x_QGP_lab_bf,vector<double> &p_QGP_lab_bf,vector<double> &x_QGP_lab_af,vector<double> &p_QGP_lab_af,vector<vector<vector<double> > > & all_jet_x_1ev,vector<vector<vector<double> > > &all_jet_p_1ev);

//*ーーーーーーーーーーーーー宣言・関数finーーーーーーーーーーーーー

//*ーーーーー　取得するデータファイル　ーーーーー
//終了時刻におけるジェットの運動量
std::ofstream fin_px_DATA_out;
std::ofstream fin_py_DATA_out;
std::ofstream fin_pz_DATA_out;
std::ofstream Etot_recoiled_QGP_1ev_out;

std::ofstream theta_af_col;
std::ofstream theta_final_ts_jet;//時間を止めた時のθ分布

int main(){


//ーーーーーーーーーーーーー取得するデータファイル一覧ーーーーーーーーーーーーーーーーーー
bool data_file_get = true;//データファイルを取得するかどうか

if(data_file_get){
    //時間発展が終了したときの運動量
    fin_px_DATA_out.open("fin_px_DATA.dat");
    fin_py_DATA_out.open("fin_py_DATA.dat");
    fin_pz_DATA_out.open("fin_pz_DATA.dat");
    //1evでjetによって差っ引かれたQGPのエネルギーの総和を記録するファイル
    Etot_recoiled_QGP_1ev_out.open("Etot_recoiled_QGP_1ev.dat");
}

//todo 今だけここを特別扱い
//θのデータファイル化
theta_af_col.open("theta_af_col.dat");
theta_final_ts_jet.open("theta_final_ts_jet");
//todo


//ーーー実験の設定で使う変数の定義ーーーーー
//double t_scale_for_L;//L(素片の1辺)に対して何倍の時間まで見るか
//double ts_max;//タイムステップループの最大値
//int resolution;//素片に対してジェットが直線に進んだ時、時間(素片の長さ)を何分割するか

//ーーー実験の設定（イベント数・最大時間・時間の解像度・seed）ーーー
cout << "events times" << endl;
//cin >> events;



//cout << "素片の1辺Lの何倍(整数値)のスケールの時間まで観察するか(大きい値をいれれば問題なし:ジェットが箱から出た時点で終了するため)" << endl;
//cin >> t_scale_for_L;//ジェットが衝突をせずに最短で辿り着いた時の時間を1とする
//*まずは観察する時間の長さを固定する

//t_scale_for_L = 1;

//cout << "素片の大きさに対して何分割するか(回分割)" << endl;
//cin >> resolution;
//*まずは解像度を固定する

//!箱から出ないように箱のサイズを決定しているので、箱の何分の一とか考えずに直接手で与える。
//!一番初めにグローバル変数として諸々定義済み

//*resolution = 10.0;
//*ts_wid=L/resolution;//ts_widは関数でも使うので既にグローバル変数として定義済み

//cout << "タイムステップ幅(fm/c=fm):"<< ts_wid * h_barc << endl;
//cout << "流体素片の1辺の長さ(fm):"<< L*h_barc << endl;

//ts_max = resolution*t_scale_for_L;//タイムステップの最大値
//cout << "最大タイムステップ数:" << ts_max << endl;

//todo A:New:jetの初期値のインプット

//*ーーーーーーーーーーーーー乱数系ーーーーーーーーーーーーー
//分布関数の範囲に積分関数を使用しているのでmain関数の中で定義する必要がある。
p_dist_bose.param(std::uniform_real_distribution<>::param_type(F_integ(xmin,0,T),F_integ(xmax,0,T)));
//std::uniform_real_distribution<double> p_dist_fermi(F_integ(xmin,1),F_integ(xmax,1));//運動量P a≤f(P)≤b
p_dist_fermi.param(std::uniform_real_distribution<>::param_type(F_integ(xmin,1,T),F_integ(xmax,1,T)));
//?　これって積分する必要あるんだっけ？普通に解析会でよくない？それならどこでも宣言できる
//→だけどいずれmassがあることを考えるなら結局積分しないといけなくなる。

//*ーーーーーーーーーーーーー全実験を通して使う配列ーーーーーーーーーーーーー
//eventsやts_maxはmain関数の中で決まる値なのでこのタイミングで宣言する
vector<int> col_total_vec_in_1ev(events);//実験ごとの衝突回数を保存

cout << "ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー結果ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー" << endl;

//eventのループ開始ーーーーーーーーーーーーーーーーーーーー
for(int t=0; t< events;t++){
    
    cout << "--------" << t+1 << "events目--------" << endl;//0+1=1イベント目からスタート

    //*ーーーーーーーーイベントごとに初期化する配列ーーーーーーーー
    col_times_total_1ev =0;//t回目の実験での衝突回数
    num_of_jet=1;//Newfin:jetの個数    //*始めはジェットの数を1つに設定する
    //tの次元,IDの次元,座標or運動量の次元で3次元配列。
    vector<vector<vector<double> > > all_jet_x_1ev(ts_max+2,vector<vector<double> >(num_of_jet,vector<double>(4)));
    vector<vector<vector<double> > > all_jet_p_1ev(ts_max+2,vector<vector<double> >(num_of_jet,vector<double>(4)));
    
    //*ーーイベント毎のデータ書き込みや保存量確認に使用する配列や変数ーー
    num_of_col_ts_1ev=0;//衝突の起きるtsの数
    num_of_no_col_ts_1ev=0;//衝突の起きないtsの数
    col_ts_ID.resize(ts_max);//衝突が起きるtsを記録しておく
    check_1ev=0;//保存量の保存を満たしているイベント数をカウント
    vector<vector<int> > inside_jet_ID(ts_max+2,vector<int>(0));//箱の内側にいるjetを記録
    vector<vector<int> > outside_jet_ID(ts_max+2,vector<int>(0));//あるts_labにおける箱の内側にいるjet_IDを記録。
    //1evにおけるQGPから差っ引かれるエネルギーの総和
    double Etot_recoiled_QGP_1ev = 0.0;

    //todo C:Newfin:複数jetの初期値を入力・取得
    //初めは手で与えて1つのjetの反跳を記述するので飛ばす。
    // for (int i = 0; i < num_of_jet ; i++){//jetの数(jetのID)のループ
    //     for (int i = 0; j < 4; j++){//座標と運動量の空間成分
    //         all_jet_x_1ev[0][i][j]=******;//[0](初期時刻t=0)の[0]IDにiを代入
    //         all_jet_p_1ev[0][i][j]=******
    //     }
    // }

    //ーーーーーーーーーーーーー流体素片に注入するジェットの初期値ーーーーーーーーーーーーー
    //時刻0(初期値)における1つ目(ID=0)のジェットの情報
    //ージェットの初期座標ー
    all_jet_x_1ev[0][0][0]= 0.0;//時間成分 fm
    all_jet_x_1ev[0][0][1]= 0.0;//x fm
    all_jet_x_1ev[0][0][2]= 0.0;//y
    all_jet_x_1ev[0][0][3]= 0.0;//z
    //ージェットの初期運動量ー
    //*(massless)
    all_jet_p_1ev[0][0][0]= ini_jet_E;//エネルギー MeV
    all_jet_p_1ev[0][0][1]= ini_jet_P;//px MeV
    all_jet_p_1ev[0][0][2]= 0.0;//py
    all_jet_p_1ev[0][0][3]= 0.0;//pz

    //タイムステップのループーーーーーーーーーーーー
    for(int ts_lab=0;ts_lab<ts_max;ts_lab++){//0→1番目のtsに進むときに衝突する可能性もあるので、ts_labは0からスタート

        cout << "ーーーーーーーts:" << ts_lab << ":ーーーーーーーーーーーーー" << endl;

        for(int i=0;i<num_of_jet;i++){//全ジェットの時間を0に初期化
            all_jet_x_1ev[ts_lab][i][0]=0.0;
            //todo C ↑リセットのやり方。宣言し直して初期化とかもできる？ただ怖いから避ける。
        }

        //タイムステップごとに初期化する配列
        coll_times_w_in_ts=0;//Newfin:jetとQGPの被りを弾く前。ここから被りのあるペアを削る。
        coll_times_in_ts=0;//Newfin:jetの被りを弾いた後。＝衝突回数
        no_col_ts=0;//衝突が起きないタイムステップかどうかの情報を入れておく。毎タイムステップ更新する。
        no_coll_jet_ID.resize(0);//衝突しないIDを初期化

        //* ーーー　割り当て　ーーー
        QGP_assign(all_QGP_x_1ts,all_QGP_p_1ts);//更新されてほしいベクトルを入れる
        //* ーーーーーーーーーーーーー全てのQGPパートンへの(x,p)の割り当て完了ーーーーーーーーーー

        if(N_tot > 0){//素片中に粒子があるなら各jetとQGPパートンの判定を行う
        //*ーーーーー　判定　ーーーーーー
            col_judge(coll_pairs_w,all_jet_x_1ev,all_jet_p_1ev,ts_lab);

            //ーー　衝突判定の結果を調べる　ーー
            if(coll_times_w_in_ts > 0){//1つでも衝突を起こすjetとQGPのペアがあるなら
                num_of_col_ts_1ev +=1;//衝突が起きるtsをカウント
                col_ts_ID[num_of_col_ts_1ev-1] = ts_lab;//衝突が起きるtsを記録しておく
                w_erase(coll_pairs_w,coll_pairs_1on1);//*被りを除去する関数
            //*被り排除パートfinーーーーー
        //*ーーーーー　判定fin　ーーーーー
            //ーーーーーーーーーーーーーーーー衝突を記述ーーーーーーーーーーーーーーーー
             //衝突が起きる回数分のループ
                for (int i = 0; i < coll_times_in_ts; i++){
                    col_times_total += 1;//全ての衝突にカウント
                    col_times_total_1ev += 1;//1イベントにおける衝突にカウント

                    //ーーーーー　衝突開始　ーーーーー
                    vector<double> x_jet_lab_bf(4),p_jet_lab_bf(4);//衝突に使う配列
                //*i番目に衝突するjetとQGPパートンの情報を取得
                    //ジェットの情報を取得, coll_pairs_1on1[i][1]→coll_ID=iの衝突ペアのjetのIDを取り出す
                    for(int j=0;j<4;j++){
                        x_jet_lab_bf[j] = all_jet_x_1ev[ts_lab][coll_pairs_1on1[i][1]][j]; //iは衝突のID,jはジェットの各成分を表す,ジェットは時間の次元をもつ4次元配列
                        p_jet_lab_bf[j] = all_jet_p_1ev[ts_lab][coll_pairs_1on1[i][1]][j];
                    }
                    //QGPの情報取得, coll_pairs_1on1[i][2]→coll_ID=iの衝突ペアのQGPのIDを取り出す
                    for(int j=0;j<4;j++){
                        x_QGP_lab_bf[j] = all_QGP_x_1ts[coll_pairs_1on1[i][2]][j]; //iはQGPのID,jはQGPの各成分を表す,QGPはts毎にサンプリングし直すので時間の次元を持たない2次元配列
                        p_QGP_lab_bf[j] = all_QGP_p_1ts[coll_pairs_1on1[i][2]][j];
                    }
                    //cout << "衝突の起こるペア j:q " << coll_pairs_1on1[i][1] << ":" << coll_pairs_1on1[i][2] << endl;
                    Etot_recoiled_QGP_1ev += p_QGP_lab_bf[0];//QGP→jetになる（jetによって差っ引かれる）分のエネルギーの総和を記録する
                    Etot_recoiled_QGP_allev += p_QGP_lab_bf[0];
                
                //*ーーーーー重心系への変換開始-----
                    //衝突前) ローレンツ変換後 エネルギーも変化するので4成分用意
                    vector<double> p_1_bf_CM(4),p_2_bf_CM(4),beta_CM_pole(3),p_1_af_CM(4),p_2_af_CM(4);
                    //衝突後の最終状態
                    vector<double>x_jet_lab_af(4),p_jet_lab_af(4),x_QGP_lab_af(4),p_QGP_lab_af(4);
                    trans_to_CM(p_jet_lab_bf,p_QGP_lab_bf,p_1_bf_CM,p_2_bf_CM,beta_CM_pole);

                    //等方散乱で角度を割り当てる
                    vector<double> p_1_bf_CM_pole(3);//座標回転のための極座標への変換
                    vector<double> p_2_bf_CM_pole(3);
                    Cart_to_Pole(p_1_bf_CM_pole,p_1_bf_CM);//直交座標を極座標に変換
                    Cart_to_Pole(p_2_bf_CM_pole,p_2_bf_CM);

                //*ーーーーーーーーーーーーー衝突ーーーーーーーーーーーーー
                    vector<double> p_1_af_CM_pole(3),p_2_af_CM_pole(3);//衝突後の極座標表示
                    //ランダムにθとφを決定する: 逆向きに割り当てるp1(|p1|,θ,φ),p2(|p2|,-θ,-φ)
                    std::uniform_real_distribution<double> dis_theta(-1,1);//-1≤cosθ≤1
                    std::uniform_real_distribution<double> dis_phi(0,2*M_PI);//0≤φ"≤"2π

                    cout << ">>>>> 衝突 <<<<<" << endl;
                    //エネルギーは変化なし:E_p1_bf_CM=E_p1_af_CM, 向きが変わるだけ:|p1_bf_CM|=|p1_af_CM|
                    //p_pole(|p|,θ,φ), 運動量の大きさ(=エネルギー(massless))は変化なし
                    p_1_af_CM_pole[0]=p_1_bf_CM_pole[0];
                    p_2_af_CM_pole[0]=p_2_bf_CM_pole[0];
                    p_1_af_CM[0]=p_1_bf_CM[0];
                    p_2_af_CM[0]=p_2_bf_CM[0];
                    //θ,φの割り当て
                    p_1_af_CM_pole[1] = acos(dis_theta(mt)) + p_1_bf_CM_pole[1];//θの割り当て: 逆向きに散乱　θ→π-θ回転
                    p_2_af_CM_pole[1] = M_PI-p_1_af_CM_pole[1] ;
                    p_1_af_CM_pole[2] = dis_phi(mt) + p_1_bf_CM_pole[2]; //φの割り当て: 逆向きに散乱　φ→φ+π回転
                    p_2_af_CM_pole[2] = p_1_af_CM_pole[2] + M_PI;

                    //直交座標に変換し直す。
                    Pole_to_Cart(p_1_af_CM,p_1_af_CM_pole);
                    Pole_to_Cart(p_2_af_CM,p_2_af_CM_pole);
                //*ーーーーーー衝突終了ーーー
                //*ーーーーーー実験室系への逆変換ーーーーーー
                    reverse_trans_to_CM(p_1_af_CM,p_2_af_CM, p_jet_lab_af, p_QGP_lab_af,beta_CM_pole);
                
                    if(conserve_check){//*保存量確認
                        conserve_check_func(p_jet_lab_bf,p_QGP_lab_bf,p_jet_lab_af,p_QGP_lab_af);
                    }

//todo θ分布の取得

/*
*p_jet_lab_bfとp_jet_lab_af, p_QGP_lab_afのなす角度を求める

acos((a・b)/(|a||b|))

space成分を移し替え

acos(in_pro_3d(pj_lab_bf_space,pj_lab_af_space)/(abs(pj_lab_bf_space)*abs(pj_lab_af_space))
*/

//space成分の移し替え
vector<double> pj_lab_bf_space(3),pj_lab_af_space(3),pQ_lab_af_space(3);

for(int i=0;i<3;i++){
    pj_lab_bf_space[i] = p_jet_lab_bf[i+1];
    pj_lab_af_space[i] = p_jet_lab_af[i+1];
    pQ_lab_af_space[i] = p_QGP_lab_af[i+1];
}

//衝突後のジェットとの角度 
//*ラジアンで表す。
theta_af_col << acos(in_pro_3d(pj_lab_bf_space,pj_lab_af_space) / ( sqrt( in_pro_3d(pj_lab_bf_space,pj_lab_bf_space) ) * sqrt( in_pro_3d(pj_lab_af_space, pj_lab_af_space) ) ) ) << endl;
//衝突後のQGPとの角度
theta_af_col << acos(in_pro_3d(pj_lab_bf_space,pQ_lab_af_space) / ( sqrt( in_pro_3d(pj_lab_bf_space,pj_lab_bf_space) ) * sqrt( in_pro_3d(pQ_lab_af_space, pQ_lab_af_space) ) ) ) << endl;



//todo θ分布の取得


                //*ーーーーー散乱後のデータ更新
                    update_af_col(i,ts_lab,x_jet_lab_bf,p_jet_lab_bf,x_jet_lab_af,p_jet_lab_af,x_QGP_lab_bf,p_QGP_lab_bf,x_QGP_lab_af,p_QGP_lab_af,all_jet_x_1ev,all_jet_p_1ev);

                    cout << "" << endl;
                    cout << "num_of_jet:" << num_of_jet << endl;
                  
                }//各衝突を記述するループfin(i→coll_times_in_ts)
                //衝突しないジェットがあるなら
                //衝突しないジェットの軌跡を進める
                cout << "no_col_jet_ID.size():" << no_coll_jet_ID.size() << endl;
                if(no_coll_jet_ID.size()>0){
       
                    for(int i=0;i<no_coll_jet_ID.size();i++){

                        //cout << "衝突しないジェットのID:" << no_coll_jet_ID[i] << endl;

                        //次のタイムステップ=ts_lab+1に値を移し替える 値の更新
                        //x_0成分はタイムステップごとに初期化するので入れる意味なし。
                        //[no_coll_jet_ID[i]]に衝突しない粒子IDが入っているので、このIDを持つものを順番に取り出していけばいい。
                        for(int j=0;j<3;j++){//空間成分を1進める
                            //cout <<"衝突しないジェットの軌跡をすすめたもの:"<< all_jet_x_1ev[ts_lab][no_coll_jet_ID[i]][j+1] + ( ts_wid * (all_jet_p_1ev[ts_lab][no_coll_jet_ID[i]][j+1]/all_jet_p_1ev[ts_lab][no_coll_jet_ID[i]][0]))<< endl;
                            all_jet_x_1ev[ts_lab+1][no_coll_jet_ID[i]][j+1] = all_jet_x_1ev[ts_lab][no_coll_jet_ID[i]][j+1] + ( ts_wid * (all_jet_p_1ev[ts_lab][no_coll_jet_ID[i]][j+1]/all_jet_p_1ev[ts_lab][no_coll_jet_ID[i]][0]));
                        }
                        for(int j=0;j<4;j++){//運動量は変化なし
                            //cout << "運動量は変化なし" << endl;
                            all_jet_p_1ev[ts_lab+1][no_coll_jet_ID[i]][j]=all_jet_p_1ev[ts_lab][no_coll_jet_ID[i]][j];
                        }
                    }
                }


            }//衝突がある場合のif
            
            //衝突粒子がないなら
            else{
                no_col_ts=1;//衝突が起きないタイムステップ→no_col_tとして=1を代入
                num_of_no_col_ts_1ev +=1;//衝突が起きないtsの数をカウント
                //cout << "ーーーーーーーーーーーーー　衝突なし　ーーーーーーーーーーーーー" << endl;
            }//衝突が起こらない場合(else)
        }//N_tot > 0

        //N_tot == 0(素片中にQGPがない)なら衝突は起こらない
        else if(N_tot==0){
            cout << "衝突なし → pを変えずに軌跡を1つ進める" << endl;
            no_col_ts=1;//衝突が起きないタイムステップ→no_col_tとして=1を代入
        }//N_tot == 0,次のタイムステップへ

        if(no_col_ts==1){//no_col_tsが1→衝突が起きないタイムステップなら何も変化させずにジェットの軌跡を進める
            //cout << "num_of_jet:" << num_of_jet << endl;
            for (int i = 0; i < num_of_jet; i++){
                //i番目のジェットの軌跡を順番に進める, 衝突がない場合は運動量を変化させずに軌跡を1つ進める, 次のタイムステップ=ts_lab+1に値を移し替える           
                for(int j=0;j<3;j++){//空間成分(j+1)を1進める, x_0成分はタイムステップごとに初期化するので入れる意味なし
                    all_jet_x_1ev[ts_lab+1][i][j+1] = all_jet_x_1ev[ts_lab][i][j+1] + ( ts_wid * (all_jet_p_1ev[ts_lab][i][j+1]/all_jet_p_1ev[ts_lab][i][0]));    
                }
                for(int j=0;j<4;j++){//運動量は変化なし
                    all_jet_p_1ev[ts_lab+1][i][j]=all_jet_p_1ev[ts_lab][i][j];
                }

            }

        }

    }//タイムステップループts_labの終わり。

    //*ーーーーーーーーーーーーーそのtsが終わる瞬間にそれぞれのジェットの位置を確認ーーーーーーーーーーーーー

    for(int j=0;j<num_of_jet;j++){

        if(data_file_get){

            fin_px_DATA_out << all_jet_p_1ev[ts_max][j][1]<< endl;
            fin_py_DATA_out << all_jet_p_1ev[ts_max][j][2]<< endl;
            fin_pz_DATA_out << all_jet_p_1ev[ts_max][j][3]<< endl;                
        }

//todo A 実験終了後のθ分布

        //space成分の移し替え
        vector<double> pj_lab_fin_space(3);
        for(int i=0;i<3;i++){
            pj_lab_fin_space[i] = all_jet_p_1ev[ts_max][j][i+1];
        }

        //衝突後のジェットとの角度 
        //p_fin_x/|p_fin|
        //*ラジアンで表す。
        theta_final_ts_jet << acos((pj_lab_fin_space[0]) / (sqrt( in_pro_3d(pj_lab_fin_space, pj_lab_fin_space) ) ) ) << endl;

        if( all_jet_x_1ev[ts_max][j][1] < -(L/2.0) || L/2.0 <all_jet_x_1ev[ts_max][j][1] ||
            all_jet_x_1ev[ts_max][j][2] < -(L/2.0) || L/2.0 <all_jet_x_1ev[ts_max][j][2] ||
            all_jet_x_1ev[ts_max][j][3] < -(L/2.0) || L/2.0 <all_jet_x_1ev[ts_max][j][3]){//箱の外に出てしまっているのなら
            num_of_out_jet_fin += 1;//素片の外にでたジェットとしてカウント
        }
    }
    cout << "終了時刻になったので実験終了" << endl;
    cout << "num_of_jet:" << num_of_jet << endl;         
    num_of_jet_fin = num_of_jet;

//*ーーーイベント終了後の諸計算ーーー
    if(check_1ev > 0){
        check_all_ev += 1;//イベントの最後に保存量をカウント
    }
    if(data_file_get){
    //このイベントで差っ引かれたQGPのエネルギーの総和
        Etot_recoiled_QGP_1ev_out << Etot_recoiled_QGP_1ev << endl;
    }
    cout << "衝突が起きるts:" << endl;
    for(int i=0;i<num_of_col_ts_1ev;i++){
        cout << "ts_ID:" << col_ts_ID[i] << endl;
    }
}//イベントループeventsの終わり

//*ーーーーーーーーーーーーー実験データの集計・表示ーーーーーーーーーーーー
cout << "" << endl;
cout << "保存量の確認" << endl;
cout << "in all cols: " << check_all_col << "/" << col_times_total << endl;

cout << "素片の外に出た粒子の数:" << num_of_out_jet_fin << endl;

cout << "" << endl;
cout << "###   実験の情報   ###" << endl;
cout << "素片の1辺の長さ(fm):" << L * h_barc << endl;
cout << "平均粒子数(個):" << N_total_mean << endl;
cout << "平均粒子数密度(個/fm^3):" << N_total_mean/pow(L*h_barc,3.0) << endl;
cout << "" << endl;

cout << "QGPのエネルギー(MeV):"<< E_mean << endl;
cout << "QGPのエネルギー密度(Mev)):"<< E_mean/pow(L*h_barc,3.0) << endl;
cout << "1evで差っ引かれたエネルギーの平均(MeV):" << Etot_recoiled_QGP_allev /events << endl;
cout << "" << endl;

double E_mean_1QGP = E_mean/N_total_mean;
double col_times_mean_1ev = col_times_total/events;//1ev当たりの平均衝突回数
cout << "概算パート" << endl;
cout << "①1個あたりのエネルギー(MeV):" <<  E_mean_1QGP << endl;
cout << "②1evあたりの衝突回数の合計(回):" << col_times_mean_1ev << endl;
cout << "①*②=1evあたりの差っ引かれたエネルギーの概算値(MeV):" << E_mean_1QGP*col_times_mean_1ev << endl;

cout << "この実験におけるジェットの初期値:" << ini_jet_E/1000.0<< " GeV" << endl;
cout << "ts_max:" << ts_max << endl;


cout << "イベント数:" << events << endl;
cout << "衝突回数:" << col_times_total << endl;

}//main関数の終わり






//ーーーーーーーーーーーーーーーーーーーーー//ーーーーーーーーーーーーーーーーーーーーー
//*ーーーーーー　　　　　　　　　main fin           ーーーーーーーーーーーーーーー
//ーーーーーーーーーーーーーーーーーーーーー//ーーーーーーーーーーーーーーーーーーーーー






//割り当てる関数
void QGP_assign(vector<vector<double> > & QGP_x,vector<vector<double> > & QGP_p){
//!includeファイルにするのは厳しそう。
//→フェルミ分布とかやらないといけないし、そのためには温度も自由度も全て引数として渡さなければいけないので。
//todo B もし作るなら
//引数 N_tot,T,d
//でも基本的に乱数は他の機構にも使うと思うから割り当てだけで完結するのは難しそう。
//おそらくプログラムは完全な状態じゃないと独立したファイルにできないので。

//二分法に必要な変数
double high_defo = xmax ;//初期値(右端)
double low_defo = 0.0;//初期値(左端)　
double cri = pow(0.1,5.0);//判定基準範囲（二分法）
double y_k;//F(P)の乱数

int num_par,num_bose,num_fermi;

    vector<double> part_data(9);//割り当てに使用。{x,y,z,E,|p_x|,|p_y|,|p_z|,i,type}; i=ID 0 <= i <= N(粒子数)
    //各QGPパートンへの(x,p)の割り当て開始ーーーーー
    //粒子タイプのループ(gluon=0,quark=1)
    for(int type=0;type<num_type_max;type++){

        //粒子数を決める
        if(type==0){
            num_par = dist_bose(mt);//bose distribution
            num_bose = num_par;
            N_tot = num_bose;

        }
        
        else if(type==1){
            num_par = dist_fermi(mt);//fermi distribution
            num_fermi = num_par;
            N_tot = num_bose + num_fermi;//2周目はboseとfermiの合計

        } 

        if(num_par>0){
            //ーーーーーーーーーーーーー割り当てpartーーーーーーーーーーーーー
            //1tsにおける全てのQGPのxとpを入れとく配列を用意
            //粒子数(N_tot)分用意する
            all_QGP_x_1ts.resize(N_tot,vector<double>(4));
            all_QGP_p_1ts.resize(N_tot,vector<double>(4));
            //粒子数分ループで座標と運動量をそれぞれのQGPに割り当てる
            for(int n=0;n<num_par;n++){

                //パラメータ(p,x)割り当て

                //時間成分は0を代入(タイムステップが切り替わった瞬間を0とおく)
                x_QGP_lab_bf[0]=0.0;
                //ーーーーー座標の割り当てーーーーー
                for (int k = 0; k < 3; ++k) { 

                    //割り当て用の配列は共用
                    part_data[k] = distr(mt);
                    //x_QGP_lab_bf[k+1]= part_data[k];//座標の空間成分

                    //全ての粒子の情報を含んだ配列ではtype別に分ける
                    if(type==0){//1周目
                        //all_x_in_1evは3成分。t_ste_labという時間成分をすでに持っているので。
                        //all_x_in_1ev[ts_lab][n][k] = part_data[k];
                        //all_QGP_x_1ts(time step)は4成分                
                        all_QGP_x_1ts[n][k+1] = part_data[k];//空間成分
                    
                    }

                    else if(type==1){//2周目
                        //all_x_in_1ev[ts_lab][num_bose + n][k] = part_data[k]; 
                        //すでにnum_bose分QGPが入っている           
                        all_QGP_x_1ts[num_bose + n][k+1] = part_data[k];
                    }
                }//ーーーーー座標の割り当てfinーーーーー
                
                //二分法による運動量の割り当て
                //初期化
                double center;
                double gap;
                double high = high_defo;
                double low = low_defo;
                //運動量の大きさの決定
                if(type==0){
                    y_k = p_dist_bose(mt);
                }
                else if(type==1){
                    y_k = p_dist_fermi(mt);   
                }

                for(int i=0;i<30;i++){ //二分法で乱数(y_k)に対応する解Pを求める
                    gap= abs(high-low);
                    center = (high + low)/2;
                    if (func (center,type,y_k,T)>=0){
                        high = center;
                    }
                    if (func (center,type,y_k,T)<0){
                        low = center;
                    }
                    if(func (center,type,y_k,T)==0||cri >= gap){

                        break;

                    }
                }
                //二分法の割り当てfin
                //r,θ,φの割り当て                        
                double r = center;// F(p)から二分法によってもとまったp
                double theta = acos(dis_theta(mt));
                double phi = dis_phi(mt) ;
                //*{x,y,z,E,|p_x|,|p_y|,|p_z|,i,type}
                //p_x,p_y,p_zの割り当て
                part_data[4] = r *sin(theta)*cos(phi);//|p_x|=rsinθcosφ
                part_data[5] = r *sin(theta)*sin(phi);//|p_y|=rsinθsinφ
                part_data[6] = r *cos(theta);//|p_z|=rcosθ
                //運動量をもとにmassless粒子のエネルギーを入れる
                part_data[3] = sqrt(pow(part_data[4],2.0)+pow(part_data[5],2.0)+pow(part_data[6],2.0));
                //違う型に代入するときは(型) 変数名 でキャストする
                part_data[8] = (double)type;//粒子の種類(gluon0,quark1)

                //割りふったデータをp_2に移し替え
                p_QGP_lab_bf[0]=part_data[3];//エネルギー

                for(int i=0;i<3;i++){
                    p_QGP_lab_bf[i+1]= part_data[i+4];//空間成分のみ                              
                }


                //x,y,z,"px,py,pz"←ここに入れる
                if(type==0){//1周目
                    part_data[7] = (double)n;//粒子ID
                    for(int k=0;k<4;k++){
                        all_QGP_p_1ts[n][k]=part_data[k+3];
                    }          
                }

                else if(type==1){//2周目
                    part_data[7] = (double)(n+ num_bose);//粒子ID
                    for(int k=0;k<4;k++){
                        all_QGP_p_1ts[num_bose + n][k]=part_data[k+3];                   
                    }
                }
            }//粒子数ループfin
            //*ーーーーーーー1QGPへの割り当てpart finーーーーーーーー
        }

        else{//num_par=0なら割り当てを行わない。→何もしない
            //*type=0ならtype=1のループへ。type=1なら割り当てのループから抜ける。

        }
    }//*ーーーーー粒子タイプループfinーーーーー
    return;
}


void col_judge(vector<vector<double> >& col_pair_W, vector<vector<vector<double> > > & all_jx_1ev,vector<vector<vector<double> > > & all_jp_1ev,
    int ts_lb){

double b_rel_sqr;//ローレンツ不変な衝突径数
double t_c_CM;//重心系での2粒子の最接近時間

//ーー　各jetに対して判定を行う　ーー
for(int j = 0; j <num_of_jet; j++){

    //ー　判定用のjet配列を準備　ー
    //判定用なのでjetのループごとに初期化
    //1つのtsとこのjetの判定にのみ使うので、tsとjetの情報は必要なし
    vector<double> jet_lab_judge_x(4),jet_lab_judge_p(4);//Newfin:判定用のjetの配列
    //1evにおけるjetの配列からts_labにおけるj番目のjetをとってくる
    for (int k = 0; k < 4; k++){
        
        jet_lab_judge_x[k] = all_jx_1ev[ts_lb][j][k];
        jet_lab_judge_p[k] = all_jp_1ev[ts_lb][j][k];
                        
    }
    
    //ーーーーーーー　各QGPパートンに対して判定を行う　ーーーーーーー
    for(int q = 0; q < N_tot;q++){

        //ー　判定用のQGP配列を準備　ー
        //判定用なのでQGPのループごとに初期化
        vector<double> QGP_lab_judge_x(4),QGP_lab_judge_p(4);//Newfin:判定用のjetの配列
        //1tsにおけるQGPの配列から、q番目のQGPをとってくる
        for (int i = 0; i < 4; i++){
            
            QGP_lab_judge_x[i] = all_QGP_x_1ts[q][i];
            QGP_lab_judge_p[i] = all_QGP_p_1ts[q][i];
            
        }

        //ーーーーーーーーーーーーー判定用にx,P,p,qを作るーーーーーーーーーーーーー
        //まずは衝突径数b_relを算出してb_relによる判定から行う
        vector<double> x_lab_bf(4),P_lab_bf(4),p_lab_bf(4),q_lab_bf(4);

        for(int i=0;i<4;i++){
            //相対座標
            x_lab_bf[i] = jet_lab_judge_x[i] - QGP_lab_judge_x[i];
            //P=p1+p2
            P_lab_bf[i] = jet_lab_judge_p[i] + QGP_lab_judge_p[i];
            //相対運動量p=p1-p2
            p_lab_bf[i] = jet_lab_judge_p[i] - QGP_lab_judge_p[i];
        }

        //内積の結果を新たな変数で置き換え
        double x_x,x_p,x_P,p_p,P_P,p_P,q_q,x_q;
        x_x=in_pro_4d(x_lab_bf,x_lab_bf);
        x_p=in_pro_4d(x_lab_bf,p_lab_bf);
        x_P=in_pro_4d(x_lab_bf,P_lab_bf);    
        p_P=in_pro_4d(p_lab_bf,P_lab_bf);
        P_P=in_pro_4d(P_lab_bf,P_lab_bf);
        p_p=in_pro_4d(p_lab_bf,p_lab_bf);

        //q,その他の計算結果の割り当て
        for(int i=0;i<4;i++){
            q_lab_bf[i]=p_lab_bf[i]-((p_P/P_P)*P_lab_bf[i]);
        }
        
        q_q=in_pro_4d(q_lab_bf,q_lab_bf);
        x_q=in_pro_4d(x_lab_bf,q_lab_bf);
        b_rel_sqr=0;//値の初期化
        b_rel_sqr= - x_x + (pow((x_P),2.0) / (P_P)) + (pow((x_q),2.0) / (q_q));

        //ーーーーー衝突"径数"b_relによる判定ーーーーーーー
        if(b_rel_sqr<=d_rel_sqr){

            //ーーーーーークリアしたら今後は衝突時刻tcによる判定ーーーーーーー
            //ーーーーーtcを求めるーーーーー
            double p1_p1,p2_p2,p1_p2,x_p1,x_p2,p2_p,p1_p,p1_P,p2_P;
            double t_c_CM,t_1_CM,tc_t1_CM,tc_t2_CM;

            p1_p1=in_pro_4d(jet_lab_judge_p,jet_lab_judge_p); 
            p2_p2=in_pro_4d(QGP_lab_judge_p,QGP_lab_judge_p);
            p1_p2=in_pro_4d(jet_lab_judge_p,QGP_lab_judge_p);
            x_p1=in_pro_4d(x_lab_bf,jet_lab_judge_p);
            x_p2=in_pro_4d(x_lab_bf,QGP_lab_judge_p);
            p1_p=in_pro_4d(jet_lab_judge_p,p_lab_bf);   
            p2_p=in_pro_4d(QGP_lab_judge_p,p_lab_bf);
            p1_P=in_pro_4d(jet_lab_judge_p,P_lab_bf);   
            p2_P=in_pro_4d(QGP_lab_judge_p,P_lab_bf);

            tc_t1_CM = (p1_P / sqrt(P_P)) * ((x_P*p2_p) - (x_p*p2_P))/((p1_p*p2_P)-(p2_p*p1_P));
            tc_t2_CM = (p2_P / sqrt(P_P)) * ((x_P*p1_p) - (x_p*p1_P))/((p1_p*p2_P)-(p2_p*p1_P));

            vector<double> jet_lab_judge_x_tc(4),QGP_lab_judge_x_tc(4);
            for(int i =0;i<4;i++){
                jet_lab_judge_x_tc[i] = jet_lab_judge_x[i] + ((tc_t1_CM / (p1_P / sqrt(P_P))) * jet_lab_judge_p[i]);
                QGP_lab_judge_x_tc[i] = QGP_lab_judge_x[i] + ((tc_t2_CM / (p2_P / sqrt(P_P))) * QGP_lab_judge_p[i]);
            }

            double tc_t1_lab = jet_lab_judge_x_tc[0];//jetの最接近時刻
            double tc_t2_lab = QGP_lab_judge_x_tc[0];//QGP

            //それぞれの粒子の最接近時間を2で割って、平均値を求める
            double tc_order_lab = (tc_t1_lab +tc_t2_lab) /2.0;
            vector<double> b_rel_vec(4);
            for(int i =0;i<4;i++){
                b_rel_vec[i] = x_lab_bf[i] + ((tc_t1_CM / (p1_P / sqrt(P_P))) * jet_lab_judge_p[i]) - ((tc_t2_CM / (p2_P / sqrt(P_P))) * QGP_lab_judge_p[i]);
            }
            b_rel_sqr = in_pro_4d(b_rel_vec,b_rel_vec);

            //ーーーーー衝突"時刻"tcによる判定ーーーーー
            if((0<= tc_order_lab) && (tc_order_lab<ts_wid)){
                //cout << "衝突判定クリア" << endl;
                coll_times_w_in_ts += 1 ;//jetとQGPの被りありの衝突回数をカウント

                //新しい衝突ペアを記録するためにresizeで要素数を増やす。
                coll_pairs_w.resize(coll_times_w_in_ts,vector<double>(4));
                //要素は0からなので、衝突ペアの要素番目=衝突回数-1
                //3番目の衝突は、2番目(衝突回数3-1)番目の要素に入れておく
                //衝突の情報を記録
                coll_pairs_w[coll_times_w_in_ts-1][0] = tc_order_lab;//時刻tc
                coll_pairs_w[coll_times_w_in_ts-1][1] = (double)j;//jetのID
                coll_pairs_w[coll_times_w_in_ts-1][2] = (double)q;//QGPのID
                coll_pairs_w[coll_times_w_in_ts-1][3] = coll_times_w_in_ts-1;//衝突のID。要素番目の衝突の意味。
            }

            //衝突時刻の審査落ち
            else{
                //cout << "b_relの審査落ち" << endl;
                //衝突しないQGPに関してはいらない情報なので特に何もしなくて良いはず
                //前回のプログラムのような配列にわざわざ代入する手間がなくなったので。
                //必要ないならこのelseすらいらないかも。いらない場合は削除する。

            }
        }

        //衝突径数の審査落ち
        else{
            //cout << "tcの審査落ち" << endl;
            //特に何もしなくていい

        }
    }//各パートンへの判定ループfin
}//各jetへの判定ループfin＝各ジェットへの判定完了
    return;
}

void w_erase(vector<vector<double> >& coll_pairs_w,vector<vector<double> >& coll_pairs_1on1){

    //被りを排除するパート
    //確認のため衝突ペアを全て表示する
    //cout << "並びかえ前の行列" << endl;
    //cout << "tc j Q coll_ID" << endl;
    // for(int i=0;i<coll_pairs_w.size();i++){

    //     for (int j = 0; j < 4; ++j) {//tc,j,Q,coll_IDの4成分
    //         cout << coll_pairs_w[i][j] << ' ' ;
    //     }
    //     cout << endl;

    // }
    //ーーーーー早い順に並べて、ダブりを解消するパートーーーーー
    //tcの値=coll_pairs_w[i][0]が小さい順に並べ替える
    //sortを使って並べ替える。[i][0]成分が小さい順にiが並び変わるはず。
    sort(coll_pairs_w.begin(), coll_pairs_w.end());

    // cout << "並びかえ後の行列" << endl;
    // cout << "tc j Q coll_ID" << endl;
    // for(int i=0;i<coll_pairs_w.size();i++){

    //     for (int j = 0; j < 4; ++j) {
    //         cout << coll_pairs_w[i][j] << ' ' ;
    //     }
    //     cout << endl;

    // }
    //並べ替えOK

    //並べ替え完了済みーーーーー
    //tcの小さいものから衝突ペアとして決定していく(カップル成立)
    for(int i=0;i<coll_times_w_in_ts;i++){//全てのペアに対して被りを確認

        if(i==0){//tcが一番小さいペアは確実にカップル成立
            coll_jet_ID.resize(0);//push_backだと上書きされないので、resizeで0にしておく(初期化)
            coll_times_in_ts += 1;//衝突回数をカウント
            coll_jet_ID.push_back(coll_pairs_w[i][1]);//iでも0でもどっちでもいい。衝突する粒子として追加

            //新しい衝突ペアを記録するためにresizeで要素数を増やす。
            coll_pairs_1on1.resize(1,vector<double>(4));

            //一番最初の要素に入れておく。
            coll_pairs_1on1[0][0] = coll_pairs_w[0][0];//並べ替えられた後
            coll_pairs_1on1[0][1] = coll_pairs_w[0][1];//jetのID
            coll_pairs_1on1[0][2] = coll_pairs_w[0][2];//QGPのID
            coll_pairs_1on1[0][3] = coll_times_in_ts-1;//coll_ID

        }

        //i=>2のループ
        else{//二番目以降のペアに関しては、パートナーが先に奪われていないかチェックする必要がある。

            bool wcount_flag = 0;//被りがあるなら1。ないなら0
            //どちか一方でも被りがあるなら
            //すでに決まっている衝突回数分ループ
            //衝突が確定した配列の中を探るループ
            for(int j=0;j<coll_times_in_ts;j++){

                //被りがないかを確認
                if(coll_pairs_1on1[j][1]==coll_pairs_w[i][1] ||coll_pairs_1on1[j][2]==coll_pairs_w[i][2]){

                    //no_coll_jet_ID.push_back(coll_pairs_w[i][1]);//すでに成立しているカップルとの被りがあったら、そのジェットは衝突しないものとする。
                    //親のループの添字はiなので、i番目の衝突ペアであることに気を付ける。
                    wcount_flag = 1;//被りがあったので1を代入
                    break;//すでに被りがあって衝突しないと分かったので、この時点でループを抜ける
                }

            }

            //被りが一つもなかったなら
            if(wcount_flag == 0){

                coll_times_in_ts += 1;//衝突回数をカウント
                coll_jet_ID.push_back(coll_pairs_w[i][1]);//被りが1つもないなら衝突ジェット

                //新しい衝突ペアを記録するためにresizeで要素数を増やす。
                coll_pairs_1on1.resize(coll_times_in_ts,vector<double>(4));

                //要素は0からなので、衝突ペアの要素番目=衝突回数-1
                //3番目の衝突は、2番目(衝突回数3-1)番目の要素に入れておく
                //衝突の情報を記録
                coll_pairs_1on1[coll_times_in_ts - 1][0] = coll_pairs_w[i][0];//並べ替えられているとすればi番目の要素を入れる
                coll_pairs_1on1[coll_times_in_ts - 1][1] = coll_pairs_w[i][1];//jetのID
                coll_pairs_1on1[coll_times_in_ts - 1][2] = coll_pairs_w[i][2];//QGPのID
                coll_pairs_1on1[coll_times_in_ts - 1][3] = coll_times_in_ts - 1;//衝突候補ペアのID

            }
        }
    }//*全衝突ペアに対しての被り検証ループfin
    //*ーー衝突する粒子が決定→coll_jet_IDが決定ーー

    //衝突しないジェットを振り分ける。
    //前ジェットから衝突するjetを引いたもの＝衝突しないjet
    //no_coll_jet_IDを決定する
    bool no_col_flag = 0;
    //全てのジェットに対して確かめる。
    //全ジェットの中でi番目
    for(int i=0;i< num_of_jet ;i++){

        bool no_col_flag = 0;//ループごとに宣言し直す

        //衝突する粒子と1つずつ照合
        for(int j=0;j<coll_times_in_ts;j++){

            //被りがないかを確認
            if(coll_jet_ID[j] == i){//i番目に割り当てられたジェットが衝突するジェットのj番目と一致するなら

                //すでに衝突粒子とみなされている
                no_col_flag = 1;//被りがあったので1を代入
                break;//すでに被りがあって衝突しないと分かったので、この時点でループを抜ける
            }

        }

        //衝突粒子の中に被りが一つもなかったなら
        if(no_col_flag == 0){

            no_coll_jet_ID.push_back(i);//衝突"しない"粒子として記録
        
        }
    }

    //最後にcoll_jet_IDとno_coll_jet_IDの数を足して確認
    cout << "size of coll_jet_ID =" << coll_jet_ID.size() << endl;
    cout << "size of no_coll_jet_ID =" << no_coll_jet_ID.size() << endl;

    //衝突ジェットのID一覧
    cout << "ーー衝突するjetのIDーー(一致するはず):" << endl;
    for(int i=0;i<coll_jet_ID.size();i++){
        //cout << coll_jet_ID[i] << endl;
        //cout << coll_pairs_1on1[i][1] << endl;
    }
    cout << "" << endl;

    //衝突しないジェットのID一覧 *衝突しない粒子があるなら表示
    if(no_coll_jet_ID.size()>1){
        //cout << "衝突しないjetのID:" << endl;
        for(int i=0;i<no_coll_jet_ID.size();i++){
            //cout << no_coll_jet_ID[i] << endl;
        }
    }
    return;
}

void update_af_col(int i,int ts_lab,vector<double> &x_jet_lab_bf,vector<double> &p_jet_lab_bf,vector<double> &x_jet_lab_af,vector<double> &p_jet_lab_af,
vector<double> &x_QGP_lab_bf,vector<double> &p_QGP_lab_bf,vector<double> &x_QGP_lab_af,vector<double> &p_QGP_lab_af,vector<vector<vector<double> > > & all_jet_x_1ev,vector<vector<vector<double> > > &all_jet_p_1ev){

    //*----散乱後のp_1とx_1を更新(元々jet)-----
    //x_1の更新
    //衝突点を求める→終値を求める→次の初期値に代入の3つのプロセス
    vector<double> x_1_tc(4);//jet側の衝突点の座標を表す配列
    //初期値→衝突点 （bfを使う)
    for(int j=0;j<3;j++){//空間成分のみ代入
        //coll_pairs_1on1[i][0]にこの衝突ペアの再接近時間が入っている, つまりtc_order_labが入っている
        x_1_tc[j+1] = x_jet_lab_bf[j+1] + ( coll_pairs_1on1[i][0] * (p_jet_lab_bf[j+1]/p_jet_lab_bf[0]));
    }
    //衝突点→終値　（afを使う）
    x_jet_lab_af[0]=ts_wid;//ts_widが経過した時がそのtsにおける時刻の終値
    for(int j=0;j<3;j++){//空間成分のみ代入
        //衝突点x_1_tc→終値af
        x_jet_lab_af[j+1] = x_1_tc[j+1] + ( (ts_wid- coll_pairs_1on1[i][0]) * (p_jet_lab_af[j+1]/p_jet_lab_af[0]));
    }

    //*-----散乱後のp_2とx_2を更新(元々QGP)-----
    //x_2の更新
    vector<double> x_2_tc(4);//jet側の衝突点の座標を表す配列
    //初期値→衝突点 （bfを使う)
    for(int j=0;j<3;j++){//空間成分のみ代入
        //初期値bf→衝突点x_2_tc
        x_2_tc[j+1] = x_QGP_lab_bf[j+1] + ( coll_pairs_1on1[i][0] * (p_QGP_lab_bf[j+1]/p_QGP_lab_bf[0]));
    }
    //衝突点→終値　（afを使う）
    x_QGP_lab_af[0]=ts_wid;
    for(int j=0;j<3;j++){//空間成分のみ代入
        //衝突点x_2_tc→終値af
        x_QGP_lab_af[j+1] = x_2_tc[j+1] + ( (ts_wid- coll_pairs_1on1[i][0]) * (p_QGP_lab_af[j+1]/p_QGP_lab_af[0]));                    
    }

    //*ーーーデータを次のtsの初期値へ代入ーーー
    //*1.元々ジェットだった粒子
    //coll_pairs_1on1[i][1]にはジェットのIDが入っている。
    for(int j=0;j<4;j++){
        //次のタイムステップの配列に追加
        all_jet_x_1ev[ts_lab+1][coll_pairs_1on1[i][1]][j]=x_jet_lab_af[j];     
        all_jet_p_1ev[ts_lab+1][coll_pairs_1on1[i][1]][j]=p_jet_lab_af[j];
    }

    //*2.QGPからジェットになった粒子　ジェットを一つ増やす
    num_of_jet +=1;//ジェットが1つ増える
    //新しくQGP→jetに追加するために増えた分だけresize, iは衝突のID, ここは0で固定する。元々時間の時限はt_maxまですでに増えているので。
    for(int rest_ts=ts_lab;rest_ts<ts_max+1;rest_ts++){//ジェットが増える分の配列を用意
        //増えたジェットの情報を入れる配列を残りのts(rest_ts)分用意。+2は予備
        all_jet_x_1ev[rest_ts].resize(num_of_jet);//新たにjetが増える分の配列を追加
        all_jet_x_1ev[rest_ts][num_of_jet-1].resize(4);     
        all_jet_p_1ev[rest_ts].resize(num_of_jet);//新たにjetが増える分の配列を追加
        all_jet_p_1ev[rest_ts][num_of_jet-1].resize(4);           
    }

    for(int j=0;j<4;j++){
        all_jet_x_1ev[ts_lab+1][num_of_jet-1][j]=x_QGP_lab_af[j];//num_of_jet-1番目のジェットとして加える。(1個目のIDは0なので)
        all_jet_p_1ev[ts_lab+1][num_of_jet-1][j]=p_QGP_lab_af[j];      
    }

}