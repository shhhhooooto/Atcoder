
class GaussLaguerre{
private:

public:
  GaussLaguerre();
  ~GaussLaguerre();

  static void Gauss38(double xini,double xfin,double *xn,double *wn);
  static void GaussLag(double xini,double xfin,double *xn,double *wn);
};
