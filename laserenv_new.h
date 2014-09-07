#ifndef LASERENV_H
#define LASERENV_H

#include "geomutils.h"

#include <complex>
#include <vector>
#include <QWidget>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
using namespace std;

typedef complex<double> cc;
typedef vector< vector<cc> > table;

class LaserEnv : public QWidget
{
   Q_OBJECT
   const int n;
   int circleNumber;
   int sz;
   int w;
   int ray1X;
   int ray1Y;
   int ray2X;
   int ray2Y;
   double deviation;
   int averageIterations;
   int motionSpeed;
   int frameX;
   int frameY;
   int frameWidth;
   int frameHeight;
   int pointY;
   QLabel* lblRay;
   QLabel* lblRayText;
   QLabel* lblLattice;
   QLabel* lblLatticeText;
   QLabel* lblRay1XY;
   QLabel* lblRay2XY;
   QLabel* lblDiffraction;
   QLabel* lblDiffractionText;
   QLabel* lblFFT;
   QLabel* lblFFTText;
   QLabel* lblSize;
   QLabel* lblWidth;
   //QLabel* lblDeviation;
   QLabel* lblIterations;
   QLabel* lblMotion;
   //QLabel* lblFrameX;
   //QLabel* lblFrameY;
   //QLabel* lblFrameHeight;
   //QLabel* lblFrameWidth;
   QLabel* lblPointY;
   QLineEdit* editSize;
   QLineEdit* editWidth;
   QLineEdit* editRay1X;
   QLineEdit* editRay1Y;
   QLineEdit* editRay2X;
   QLineEdit* editRay2Y;
   //QLineEdit* editDeviation;
   QLineEdit* editIterations;
   QLineEdit* editMotion;
   //QLineEdit* editFrameX;
   //QLineEdit* editFrameY;
   //QLineEdit* editFrameHeight;
   //QLineEdit* editFrameWidth;
   QLineEdit* editPointY;
   QPushButton* iterateButton;
   QPushButton* generateButton;
   QPushButton* motionButton;
private:
   void  fftrec(vector<cc>& a, bool inv);
   table fourierTransform2D(const table& t, bool invert);
   table fastFourierTransform2D(const table& t, bool invert);
   void  print(const table& t);
   table generateBaseRay(int cx, int cy);
   QImage drawCircularLattice(int canvasWidth);
   void processCircleCenter(QImage& img, int I, int J, vector<int>& circ);
   table multiplyTables(const table& t1, const table& t2);
   table read(const QImage& img);
   QImage write(const table& t);
   QImage writeSpectrumImage(const vector<cc>& row);
   void createGUIControls();
   void setupGUIControls();
   QLabel* createLabel(const string& text);
   double calculateRayPhi(double, double, double, double);
public:
   LaserEnv(QWidget* = 0);
private slots:
   void iterateClicked();
   void generateClicked();
   void motionClicked();
};

#endif // LASERENV_H
