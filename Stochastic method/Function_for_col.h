#ifndef FUNCTION_FOR_COL
#define FUNCTION_FOR_COL

//関数内で完結するやつのみ移行した
//粒子すうやtsによって変わるやつなどは引数が増えすぎてしまうので、シンプルかつ汎用性の高いもののみ移行

//二分法で使用する関数
double F_integ(double p,int part_type_for_F_integ,double T);
double func (double x_2,int part_type,double y_k,double T);

//β(v/c)を入力して、ラピディティyを出力する関数
double beta_to_y(double beta);

//三次元の内積
double in_pro_3d(std::vector<double> vec_1, std::vector<double> vec_2);
//四次元の内積
double in_pro_4d(std::vector<double> & vec_1,std::vector<double> & vec_2);

void cal_i_1_matrix(std::vector<double>& C,std::vector<std::vector<double> >& A, std::vector<double>& B);

void new_y_CM(double y);

void new_rot_mat_phi(double phi);

void new_rot_mat_theta(double theta);

double Cart_to_Pole(std::vector<double>& P,std::vector<double>& C);

double Pole_to_Cart(std::vector<double>&C ,std::vector<double>& P);

void trans_to_CM(std::vector<double> & p_j_l_b, std::vector<double> p_Q_l_b,
    std::vector<double> & p1_b_CM,std::vector<double> & p2_b_CM,std::vector<double> & beta_CM_pole);

void reverse_trans_to_CM(std::vector<double> & p_1_af_CM, std::vector<double> p_2_af_CM,
    std::vector<double> & p_jet_lab_af,std::vector<double> & p_QGP_lab_af,std::vector<double> & beta_CM_pole);

void conserve_check_func(std::vector<double> & p_jet_lab_bf,std::vector<double> & p_QGP_lab_bf,
    std::vector<double> & p_jet_lab_af,std::vector<double> & p_QGP_lab_af);

#endif // UNCTION_FOR_COL