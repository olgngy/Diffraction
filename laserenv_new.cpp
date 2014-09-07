#include "laserenv.h"
#include "utils.h"
#include "geomutils.h"

#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <queue>
#include <QWidget>
#include <QMessageBox>
#include <QImage>
#include <QPainter>
#include <QLabel>
#include <QPixmap>
#include <QLineEdit>
#include <QPushButton>
using namespace std;

// For x in [0, n)
void LaserEnv::fftrec(vector<cc>& a, bool inv) {
   int n = a.size();
   if (n == 1) return;
   vector<cc> a0(n/2), a1(n/2);
   for (int i = 0, j = 0; i < n; i+=2, j++)
      a0[j] = a[i], a1[j] = a[i+1];
   fftrec(a0, inv);
   fftrec(a1, inv);
   double ang = 2*Pi/n * (inv ? -1 : 1);
   cc w(1), wn(cos(ang), sin(ang));
   for (int i = 0; i < n/2; i++) {
      a[i] = a0[i] + w*a1[i];
      a[i+n/2] = a0[i] - w*a1[i];
      if (inv) a[i] /= 2, a[i+n/2] /= 2;
      w *= wn;
   }
}

// For x in [-n/2, n/2), y in [-n/2, n/2)
table LaserEnv::fourierTransform2D(const table& t, bool invert) {
   table ret(n, vector<cc>(n));
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         ret[i][j] = cc(0);
         for (int i2 = 0; i2 < n; i2++)
            for (int j2 = 0; j2 < n; j2++) {
               double ang1 = 2*Pi*(i-n/2)*(i2-n/2)/n * (invert ? -1 : 1);
               double ang2 = 2*Pi*(j-n/2)*(j2-n/2)/n * (invert ? -1 : 1);
               ret[i][j] += t[i2][j2] * cc(cos(ang1), sin(ang1)) * cc(cos(ang2), sin(ang2));
            }
         ret[i][j] /= n;
      }
   return ret;
}

// For x in [-n/2, n/2), y in [-n/2, n/2)
table LaserEnv::fastFourierTransform2D(const table& t, bool invert) {
   table T = t;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         if (invert) T[i][j] *= n;
         else T[i][j] /= n; // Adjust
   if (invert) // Swap diagonal quadrants
      for (int i = 0; i < n/2; i++)
         for (int j = 0; j < n/2; j++) {
            swap(T[i][j], T[i+n/2][j+n/2]);
            swap(T[i][j+n/2], T[i+n/2][j]);
         }
   for (int i = 0; i < n; i++)
      fftrec(T[i], invert);
   for (int i = 0; i < n; i++) // Transpose
      for (int j = i+1; j < n; j++)
         swap(T[i][j], T[j][i]);
   for (int i = 0; i < n; i++)
      fftrec(T[i], invert);
   for (int i = 0; i < n; i++) // Transpose
      for (int j = i+1; j < n; j++)
         swap(T[i][j], T[j][i]);
   if (!invert) // Swap diagonal quadrants
      for (int i = 0; i < n/2; i++)
         for (int j = 0; j < n/2; j++) {
            swap(T[i][j], T[i+n/2][j+n/2]);
            swap(T[i][j+n/2], T[i+n/2][j]);
         }
   return T;
}

void LaserEnv::print(const table& t) {
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++)
         printf("(%+.3lf %+.3lf)\t\t", t[i][j].real(), t[i][j].imag());
   puts("");
}

table LaserEnv::generateBaseRay(int cx, int cy) {
   table ret(n, vector<cc>(n));
   //double U0 = 1;
   double maxi = 0;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         double val = exp(-((i-cx)*(i-cx) + (j-cy)*(j-cy))/(2.*w*w)) * 1000;
         maxi = max(maxi, val);
      }
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         double val = exp(-((i-cx)*(i-cx) + (j-cy)*(j-cy))/(2.*w*w)) * 1000;
         ret[i][j] = cc(val/maxi, 0);
      }
   return ret;
}

void LaserEnv::processCircleCenter(QImage& img, int I, int J, vector<int>& circ) {
   int R = int(sz/3. + (rand()*1./RAND_MAX)*sz*2./3 + 0.5);
   int delta = sz/3;
   for (int i = 0; i < (int)circ.size(); i+=3) {
      int I2 = circ[i];
      int J2 = circ[i+1];
      int R2 = circ[i+2];
      if ((I2-I)*(I2-I) + (J2-J)*(J2-J) <= (R+R2+delta)*(R+R2+delta))
         return;
   }
   circ.push_back(I);
   circ.push_back(J);
   circ.push_back(R);
}

QImage LaserEnv::drawCircularLattice(int canvasWidth) {
    circleNumber = 0;
   QImage ret = QImage(canvasWidth, n, QImage::Format_Mono);
   for (int i = 0; i < n; i++)
      for (int j = 0; j < canvasWidth; j++)
         ret.setPixel(j, i, 0);
   vector<int> p1(n), p2(canvasWidth), circ;
   for (int i = 0; i < n; i++)
      p1[i] = i;
   for (int i = 0; i < canvasWidth; i++)
      p2[i] = i;
   //random_shuffle(p1.begin(), p1.end());
   //random_shuffle(p2.begin(), p2.end());
   for (int i = 0; i < n; i++)
      for (int j = 0; j < canvasWidth; j++)
         processCircleCenter(ret, p1[i], p2[j], circ);
   QPainter painter(&ret);
   painter.setBrush(QBrush(Qt::white));
   for (int i = 0; i < (int)circ.size(); i+=3) {
      int I = circ[i], J = circ[i+1], R = circ[i+2];
      painter.drawEllipse(J-R, I-R, 2*R, 2*R);
      /*if (I > frameX && I < frameX + frameWidth && J > frameY && J < frameY + frameHeight)
          circleNumber++;
          */
   }
   //cout << circleNumber << endl;
   return ret;
}

table LaserEnv::multiplyTables(const table& t1, const table& t2) {
   table ret(n, vector<cc>(n));
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         ret[i][j] = t1[i][j] * t2[i][j];
   return ret;
}

double LaserEnv::calculateRayPhi(double CIRCLE_CENTER_X,
                                 double CIRCLE_CENTER_Y,
                                 double CIRCLE_RADIUS,
                                 double RAY_Y) {
    const double N1 = 1.0;
    const double N2 = 1.5;
    pt ray1 = pt(0, RAY_Y), ray2 = pt(1, RAY_Y);
    circle c = circle(CIRCLE_CENTER_X, CIRCLE_CENTER_Y, CIRCLE_RADIUS);
    line l1 = line(ray1, ray2);
    vector<pt> intersectionPoints = line_circle_intersect(l1, c);
    //if (intersectionPoints.size() < 2)
    //    throw 1;
    pt pInter = intersectionPoints.front();
    pt pCenter = pt(CIRCLE_CENTER_X, CIRCLE_CENTER_Y);
    pt dirNormal = (pCenter-pInter).normalize();
    pt dirRay = (ray2-ray1).normalize();
    double theta1 = rotationAngle(dirNormal*-1, dirRay*-1);
    double theta2 = asin(sin(theta1)*N1/N2) + Pi;
    pt dirRefracted = (rotate(pInter-dirNormal, theta2, pInter)-pInter).normalize();
    line l2 = line(pInter, pInter+dirRefracted);
    vector<pt> intersectionPoints2 = line_circle_intersect(l2, c);
    pt pInter2 = intersectionPoints2.back();
    pt dirNormal2 = (pInter2-pCenter).normalize();
    double angleTouch = rotationAngle(dirNormal2*-1, dirRefracted*-1);
    pt dirRefracted2 = (rotate(pInter2-dirNormal2, -angleTouch, pInter2)-pInter2).normalize();
    dirRefracted2 = dirRefracted2*-1;
    line l3 = line(pInter2, pInter2+dirRefracted2);
    vector<pt> intersectionPoints3 = line_circle_intersect(l3, c);
    pt pInter3 = intersectionPoints3.front();
    pt dirNormal3 = (pInter3-pCenter).normalize();
    pt dirRay3 = (pInter3-pInter2).normalize();
    double theta3 = rotationAngle(dirNormal3, dirRay3);
    double theta4 = asin(sin(theta3)*N2/N1) + Pi;
    pt dirRefracted4 = (rotate(pInter3-dirNormal3, theta4, pInter3)-pInter3).normalize();
    line l4 = line(pInter3, pInter3+dirRefracted4);
    pt pFinish;
    line_line_intersect(l4, line(pt(0, 0), pt(0, 1)), pFinish);
    double L = (pInter2-pInter).len() + (pInter3-pInter2).len();
    const double Lambda = 0.6;
    double Phi = 2*L*N2 / Lambda;
    // Adjust to [0..2Pi)
    Phi = Phi - (int)(Phi/2)*2;
    Phi *= Pi;
    //cerr << "L = " << L << endl;
    return Phi;
}

table LaserEnv::read(const QImage& img) {
   table ret(n, vector<cc>(n));
   const double Pi = acos(0.0)*2;
   vector< vector<bool> > used(n, vector<bool>(n));
   const int di[] = {-1, 0, 1, 0};
   const int dj[] = {0, 1, 0, -1};
   queue<int> q;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         QRgb rgb = img.pixel(j, i);
         bool black = !qRed(rgb) && !qGreen(rgb) && !qBlue(rgb);
         if (black || used[i][j]) continue;
         vector<int> connected;
         q.push(i);
         q.push(j);
         connected.push_back(i);
         connected.push_back(j);
         used[i][j] = true;
         int mini = 1e9, maxi = -1e9;
         int minj = 1e9, maxj = -1e9;
         while (!q.empty()) {
            int vi = q.front(); q.pop();
            int vj = q.front(); q.pop();
            for (int k = 0; k < 4; k++) {
               int ni = vi+di[k];
               int nj = vj+dj[k];
               if (ni >= 0 && ni < n && nj >= 0 && nj < n && !used[ni][nj]) {
                  QRgb rgb2 = img.pixel(nj, ni);
                  bool black2 = !qRed(rgb2) && !qGreen(rgb2) && !qBlue(rgb2);
                  if (black2) continue;
                  used[ni][nj] = true;
                  q.push(ni);
                  q.push(nj);
                  connected.push_back(ni);
                  connected.push_back(nj);
                  mini = min(mini, ni);
                  maxi = max(maxi, ni);
                  minj = min(minj, nj);
                  maxj = max(maxj, nj);
               }
            }
         }
         double r = (1+rand())*1./(RAND_MAX+1.);
         //double phi = (1+rand())*1./(RAND_MAX+1.);
         //double standardNormalDist = cos(2*Pi*phi) * sqrt(-2*log(r));
         //double normalVal = deviation * standardNormalDist;
         double cx = 0.5*(mini+maxi);
         double cy = 0.5*(minj+maxj);
         double cr = max(maxi-mini, maxj-minj)+1;
         for (int k = 0; k < (int)connected.size(); k+=2) {
            // Point inside the circle
            int vi = connected[k];
            int vj = connected[k+1];
            double distToCenterOfCircle = (pt(cx, cy)-pt(vi, vj)).len();
            double ry = cy + distToCenterOfCircle;
            //cerr << "cx = " << cx << endl;
            //cerr << "cy = " << cy << endl;
            //cerr << "cr = " << cr << endl;
            double phiVal = calculateRayPhi(cx, cy, cr, ry);
            //cerr << "phiVal = " << phiVal << endl;
            //cerr << "normalVal = " << normalVal << endl;
            ret[vi][vj] = cc(1, phiVal);
         }
      }
   return ret;
}

QImage LaserEnv::write(const table& t) {
   QImage ret = QImage(n, n, QImage::Format_Indexed8);
   QVector<QRgb> colors;
   for (int i = 0; i < 256; i++)
      colors.push_back(qRgb(i, i, i));
   ret.setColorTable(colors);
   double mini = 1e10, maxi = 0;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         double a = abs(t[i][j]);
         double v = log(a*a + 10*deviation);
         mini = min(mini, v);
         maxi = max(maxi, v);
      }
   //printf("N = %d\n", n);
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++) {
         double a = abs(t[i][j]);
         double v = log(a*a + 10*deviation);
         int col = ((v-mini)/(maxi-mini))*255;
         ret.setPixel(j, i, col);
         //printf("%d ", col);
      }
   return ret;
}

QLabel* LaserEnv::createLabel(const string& text) {
    QLabel* label = new QLabel(text.c_str(), this);
    label->setFont(QFont("Arial", 10));
    return label;
}

void LaserEnv::createGUIControls() {
    setWindowTitle("Author: Olga Nad");
    this->setFixedSize(1100, 365);
    const int d = 15;
    lblRay = new QLabel(this);
    lblRay->move(d, d);
    lblLattice = new QLabel(this);
    lblLattice->move(d+n+d, d);
    lblRay1XY = createLabel("Ray1 (x,y):");
    lblRay1XY->move(d+180, d+n+d*2+3);
    lblRay1XY->show();
    lblRay2XY = createLabel("Ray2 (x,y):");
    lblRay2XY->move(d+180, d+n+d*4+3);
    lblRay2XY->show();
    lblDiffraction = new QLabel(this);
    lblDiffraction->move(d+(n+d)*2, d);
    lblFFT = new QLabel(this);
    lblFFT->move(d+(n+d)*3, d);
    lblLatticeText = createLabel("Lattice");
    lblLatticeText->move(d+n/3+n+d, d+n+3);
    lblLatticeText->show();
    lblRayText = createLabel("Laser beam");
    lblRayText->move(d+n/3, d+n+3);
    lblRayText->show();
    lblDiffractionText = createLabel("...");
    lblDiffractionText->move(d+n/3+(n+d)*2, d+n+3);
    lblDiffractionText->show();
    lblFFTText = createLabel("Diffraction");
    lblFFTText->move(d+n/3+(n+d)*3, d+n+3);
    lblFFTText->show();
    lblSize = createLabel("Size of lattice (px):");
    lblSize->move(d, d+n+d*2+3);
    lblSize->show();
    lblWidth = createLabel("Width of beam (px):");
    lblWidth->move(d, d+n+d*4+3);
    lblWidth->show();
    /*
    lblDeviation = createLabel("Deviation:");
    lblDeviation->move(d+350, d+n+d*2+3);
    lblDeviation->show();
    */
    lblIterations = createLabel("Iterations:");
    lblIterations->move(d+830, d+n+d*4+3);
    lblIterations->show();
/*
    lblFrameX = createLabel("Frame X-coord (px):");
    lblFrameX->move(d+590, d+n+d*2+3);
    lblFrameX->show();
    lblFrameY = createLabel("Frame Y-coord (px):");
    lblFrameY->move(d+590, d+n+d*4+3);
    lblFrameY->show();
    lblFrameWidth = createLabel("Frame width (px):");
    lblFrameWidth->move(d+770, d+n+d*2+3);
    lblFrameWidth->show();
    lblFrameHeight = createLabel("Frame height (px):");
    lblFrameHeight->move(d+770, d+n+d*4+3);
    lblFrameHeight->show();
    */
    lblPointY = createLabel("Point Y-coord (px):");
    lblPointY->move(d+470, d+n+d*2+3);
    lblPointY->show();
    lblMotion = createLabel("Velocity (px):");
    lblMotion->move(d+470, d+n+d*4+3);
    lblMotion->show();
    editSize = new QLineEdit(this);
    editSize->move(d+130, d+n+d*2+1);
    editSize->resize(40, 20);
    editSize->show();
    editWidth = new QLineEdit(this);
    editWidth->move(d+130, d+n+d*4+1);
    editWidth->resize(40, 20);
    editWidth->show();
    editRay1X = new QLineEdit(this);
    editRay1X->move(d+255, d+n+d*2+1);
    editRay1X->resize(40, 20);
    editRay1X->show();
    editRay1Y = new QLineEdit(this);
    editRay1Y->move(d+300, d+n+d*2+1);
    editRay1Y->resize(40, 20);
    editRay1Y->show();
    editRay2X = new QLineEdit(this);
    editRay2X->move(d+255, d+n+d*4+1);
    editRay2X->resize(40, 20);
    editRay2X->show();
    editRay2Y = new QLineEdit(this);
    editRay2Y->move(d+300, d+n+d*4+1);
    editRay2Y->resize(40, 20);
    editRay2Y->show();
    /*
    editDeviation = new QLineEdit(this);
    editDeviation->move(d+420, d+n+d*2+1);
    editDeviation->resize(40, 20);
    editDeviation->show();
    */
    editIterations = new QLineEdit(this);
    editIterations->move(d+890, d+n+d*4+1);
    editIterations->resize(40, 20);
    editIterations->show();
/*
    editFrameX = new QLineEdit(this);
    editFrameX->move(d+720, d+n+d*2+1);
    editFrameX->resize(40, 20);
    editFrameX->show();
    editFrameY = new QLineEdit(this);
    editFrameY->move(d+720, d+n+d*4+1);
    editFrameY->resize(40, 20);
    editFrameY->show();
    editFrameHeight = new QLineEdit(this);
    editFrameHeight->move(d+890, d+n+d*4+1);
    editFrameHeight->resize(40, 20);
    editFrameHeight->show();
    editFrameWidth = new QLineEdit(this);
    editFrameWidth->move(d+890, d+n+d*2+1);
    editFrameWidth->resize(40, 20);
    editFrameWidth->show();
    */
    editPointY = new QLineEdit(this);
    editPointY->move(d+590, d+n+d*2+1);
    editPointY->resize(40, 20);
    editPointY->show();
    editMotion = new QLineEdit(this);
    editMotion->move(d+590, d+n+d*4+1);
    editMotion->resize(40, 20);
    editMotion->show();
    generateButton = new QPushButton(this);
    generateButton->setText("Generate");
    generateButton->resize(110, 25);
    generateButton->move(d+350, d+n+d*4-2);
    generateButton->show();

    iterateButton = new QPushButton(this);
    iterateButton->setText("Print average");
    iterateButton->resize(115, 25);
    iterateButton->move(d+940, d+n+d*4-2);
    iterateButton->show();

    motionButton = new QPushButton(this);
    motionButton->setText("Move");
    motionButton->resize(135, 25);
    motionButton->move(d+640, d+n+d*4-2);
    motionButton->show();
    connect(generateButton, SIGNAL(clicked()), this, SLOT(generateClicked()));
    connect(iterateButton, SIGNAL(clicked()), this, SLOT(iterateClicked()));
    connect(motionButton, SIGNAL(clicked()), this, SLOT(motionClicked()));
}

void LaserEnv::setupGUIControls() {
   // Laser beam
   table tabRay = generateBaseRay(ray1Y, ray1X);
   table tabRay2 = generateBaseRay(ray2Y, ray2X);
   table tabBoth = tabRay;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         tabBoth[i][j] = tabBoth[i][j] + tabRay2[i][j];
   QImage imgBoth = write(tabBoth);
   lblRay->setPixmap(QPixmap::fromImage(imgBoth));
   lblRay->show();
   // Circular lattice
   QImage imgLattice = drawCircularLattice(n);
   lblLattice->setPixmap(QPixmap::fromImage(imgLattice));
   lblLattice->show();
   // Diffraction
   table tabLattice = read(imgLattice);
   table tabDiffraction = multiplyTables(tabRay, tabLattice);
   table tabDiffraction2 = multiplyTables(tabRay2, tabLattice);
   table tabDiffractionBoth = tabDiffraction;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
          tabDiffractionBoth[i][j] =
                  abs(tabDiffraction[i][j])*abs(tabDiffraction[i][j]) +
                  tabDiffraction[i][j]*conj(tabDiffraction2[i][j]) +
                  conj(tabDiffraction[i][j])*tabDiffraction2[i][j] +
                  abs(tabDiffraction2[i][j])*abs(tabDiffraction2[i][j]);
   QImage imgDiffraction = write(tabDiffractionBoth);
   lblDiffraction->setPixmap(QPixmap::fromImage(imgDiffraction));
   lblDiffraction->show();
   // Diffraction's FFT
   table t4 = fastFourierTransform2D(tabDiffractionBoth, false);
   QImage imgFFT = write(t4);
   lblFFT->setPixmap(QPixmap::fromImage(imgFFT));
   lblFFT->show();
   // Write text output
   /*freopen("diffraction.txt", "wt", stdout);
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++)
         printf("%.2lf %.2lf\t", t4[i][j].real(), t4[i][j].imag());
   freopen("diffraction_abs.txt", "wt", stdout);
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++)
         printf("%.2lf, ", abs(t4[i][j]));
   */
   // Control panel
   editSize->setText(intToString(sz).c_str());
   editWidth->setText(intToString(w).c_str());
   editRay1X->setText(intToString(ray1X).c_str());
   editRay1Y->setText(intToString(ray1Y).c_str());
   editRay2X->setText(intToString(ray2X).c_str());
   editRay2Y->setText(intToString(ray2Y).c_str());
   //editDeviation->setText(doubleToString(deviation).c_str());
   editIterations->setText(intToString(averageIterations).c_str());
   /*
   editFrameX->setText(intToString(frameX).c_str());
   editFrameY->setText(intToString(frameY).c_str());
   editFrameWidth->setText(intToString(frameWidth).c_str());
   editFrameHeight->setText(intToString(frameHeight).c_str());
   */
   editPointY->setText(intToString(pointY).c_str());
   editMotion->setText(intToString(motionSpeed).c_str());
}

void LaserEnv::generateClicked() {
   sz = atoi(editSize->text().toStdString().c_str());
   w = atof(editWidth->text().toStdString().c_str());
   averageIterations = atoi(editIterations->text().toStdString().c_str());
   motionSpeed = atoi(editMotion->text().toStdString().c_str());
   ray1X = atoi(editRay1X->text().toStdString().c_str());
   ray1Y = atoi(editRay1Y->text().toStdString().c_str());
   ray2X = atoi(editRay2X->text().toStdString().c_str());
   ray2Y = atoi(editRay2Y->text().toStdString().c_str());
   /*
   frameX = atoi(editFrameX->text().toStdString().c_str());
   frameY = atoi(editFrameY->text().toStdString().c_str());
   frameHeight = atoi(editFrameHeight->text().toStdString().c_str());
   frameWidth = atoi(editFrameWidth->text().toStdString().c_str());
   */
   pointY = atoi(editPointY->text().toStdString().c_str());
   //sscanf(editDeviation->text().toStdString().c_str(), "%lf", &deviation);
   setupGUIControls();
}

void LaserEnv::iterateClicked() {
/*
   generateClicked();
   table tabRay = generateBaseRay(n/2, n/2);
   QImage imgRay = write(tabRay);
   QImage imgLattice = drawCircularLattice(n);
   table tDiffFFTAverage(n, vector<cc>(n));
   for (int i = 0; i < averageIterations; i++) {
      table tabLattice = read(imgLattice);
      table tabDiffraction = multiplyTables(tabRay, tabLattice);
      table tFFT = fastFourierTransform2D(tabDiffraction, false);
      for (int j = 0; j < n; j++)
         for (int k = 0; k < n; k++)
            tDiffFFTAverage[j][k] += cc(abs(tFFT[j][k]), 0);
      QImage imgFFT = write(tFFT);
      char fileName[64];
      sprintf(fileName, "DiffractionDFT-%d.png", i+1);
      imgFFT.save(fileName, "png");
   }
   for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
         tDiffFFTAverage[j][k] /= averageIterations;
   QImage imgDiffFFTAverage = write(tDiffFFTAverage);
   imgDiffFFTAverage.save("DiffractionDFT-AVERAGE.png", "png");
   QLabel* imageLabel = new QLabel();
   imageLabel->setPixmap(QPixmap::fromImage(imgDiffFFTAverage));
   imageLabel->setFixedSize(imgDiffFFTAverage.width(), imgDiffFFTAverage.height());
   imageLabel->show();
   // Write text output
   freopen("diffraction_average.txt", "wt", stdout);
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++)
         printf("%.2lf %.2lf\t", tDiffFFTAverage[i][j].real(), tDiffFFTAverage[i][j].imag());
   freopen("diffraction_average_abs.txt", "wt", stdout);
   for (int i = 0; i < n; i++) {
      printf("{");
      for (int j = 0; j < n; j++) {
         printf("%.2lf", abs(tDiffFFTAverage[i][j]));
         if (j != n-1) printf(", ");
      }
      printf("},\n");
   }
   //QMessageBox msgBox;
   //msgBox.setText("Resultant image generated!");
   //msgBox.exec();
*/

   generateClicked();
   table tabRay = generateBaseRay(ray1Y, ray1X);
   table tabRay2 = generateBaseRay(ray2Y, ray2X);
   table tabBoth = tabRay;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         tabBoth[i][j] = tabBoth[i][j] + tabRay2[i][j];
   QImage imgRay = write(tabBoth);
   QImage imgLattice = drawCircularLattice(n);

   cout << "n = " << n << endl;
   cout << "averageIterations = " << averageIterations << endl;
   vector<cc> motionAverage(n/averageIterations);
   for (int i = 0; i < averageIterations; i++) {
      table tabLattice = read(imgLattice);
      table tabDiffraction = multiplyTables(tabRay, tabLattice);
      table tabDiffraction2 = multiplyTables(tabRay2, tabLattice);
      table tabDiffractionBoth = tabDiffraction;
      for (int ii = 0; ii < n; ii++)
         for (int j = 0; j < n; j++)
            tabDiffractionBoth[ii][j] =
               abs(tabDiffraction[ii][j])*abs(tabDiffraction[ii][j]) +
               tabDiffraction[ii][j]*conj(tabDiffraction2[ii][j]) +
               conj(tabDiffraction[ii][j])*tabDiffraction2[ii][j] +
               abs(tabDiffraction2[ii][j])*abs(tabDiffraction2[ii][j]);
      table tabFFT = fastFourierTransform2D(tabDiffractionBoth, false);
      QImage imgFFT = write(tabFFT);
      vector<cc> meanValues;
      for (int shift = 0; shift < n; shift += motionSpeed) {
         const int i1 = pointY, j1 = shift;
         imgFFT.setPixel(j1, i1, 255);
         cc curMean;
         curMean = tabFFT[i1][j1];
         meanValues.push_back(curMean);
      }
      QImage spectrumImage = writeSpectrumImage(meanValues);
      char fileName[64];
      sprintf(fileName, "Motion-%d.png", i+1);
      spectrumImage.save(fileName, "png");
      /*QLabel* imageLabel = new QLabel();
      imageLabel->setPixmap(QPixmap::fromImage(spectrumImage));
      imageLabel->setFixedSize(spectrumImage.width(), spectrumImage.height());
      imageLabel->show();
      */
      for (int k = 0; k < n/averageIterations; k++)
         motionAverage[k] += meanValues[k];
   }
   for (int k = 0; k < n/averageIterations; k++)
      motionAverage[k] /= averageIterations;
   QImage imgMotionAverage = writeSpectrumImage(motionAverage);
   imgMotionAverage.save("!MOTION-AVERAGE.png", "png");
   QLabel* imageLabel = new QLabel();
   imageLabel->setPixmap(QPixmap::fromImage(imgMotionAverage));
   imageLabel->setFixedSize(imgMotionAverage.width(), imgMotionAverage.height());
   imageLabel->show();
}

void LaserEnv::motionClicked() {
   generateClicked();
   table tabRay = generateBaseRay(ray1Y, ray1X);
   table tabRay2 = generateBaseRay(ray2Y, ray2X);
   table tabBoth = tabRay;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         tabBoth[i][j] = tabBoth[i][j] + tabRay2[i][j];
   QImage imgRay = write(tabBoth);
   QImage imgLattice = drawCircularLattice(n);
   table tabLattice = read(imgLattice);
   table tabDiffraction = multiplyTables(tabRay, tabLattice);
   table tabDiffraction2 = multiplyTables(tabRay2, tabLattice);
   table tabDiffractionBoth = tabDiffraction;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
          tabDiffractionBoth[i][j] =
                  abs(tabDiffraction[i][j])*abs(tabDiffraction[i][j]) +
                  tabDiffraction[i][j]*conj(tabDiffraction2[i][j]) +
                  conj(tabDiffraction[i][j])*tabDiffraction2[i][j] +
                  abs(tabDiffraction2[i][j])*abs(tabDiffraction2[i][j]);
   table tabFFT = fastFourierTransform2D(tabDiffractionBoth, false);
   QImage writeLattice = write(tabDiffractionBoth);
   QImage writeFFT = write(tabFFT);
   vector<cc> meanValues;
   for (int shift = 0, it = 0; shift < n; shift += motionSpeed, it++) {
      /*
      QImage curLattice = QImage(n, n, QImage::Format_Mono);
      for (int i = 0; i < n; i++)
         for (int j = 0; j < n; j++) {
            QRgb rgb = imgLattice.pixel(j+shift, i);
            bool black = !qRed(rgb) && !qGreen(rgb) && !qBlue(rgb);
            curLattice.setPixel(j, i, !black);
         }
      */

      const int i1 = pointY, j1 = shift;
      /*
      const int i2 = frameY+frameHeight, j2 = frameX+frameWidth;
      for (int i = i1; i <= i2; i++) {
         writeFFT.setPixel(j1, i, 255);
         writeFFT.setPixel(j2, i, 255);
      }
      for (int j = j1; j <= j2; j++) {
         writeFFT.setPixel(j, i1, 255);
         writeFFT.setPixel(j, i2, 255);
      }
      */
      writeFFT.setPixel(j1, i1, 255);
      cc curMean;
      /*
      int cnt = 0;
      for (int i = i1; i <= i2; i++)
         for (int j = j1; j <= j2; j++)
            curMean += tabFFT[i][j], cnt++;
      curMean /= cnt;
      */
      curMean = tabFFT[i1][j1];
      meanValues.push_back(curMean);
      /*char fileName[64];
      sprintf(fileName, "Motion-%d.png", it+1);
      writeLattice.save(fileName, "png");
      sprintf(fileName, "Motion-Diffraction-%d.png", it+1);
      writeFFT.save(fileName, "png");
      */
   }
   /*
   freopen("motion_DFT.txt", "wt", stdout);
   printf("Initial values: ");
   double MV = 0;
   for (int i = 0; i < (int)meanValues.size(); i++) {
      double meanAbs = abs(meanValues[i]);
      MV += meanAbs;
      printf("%.12lf\t", meanAbs);
   }
   puts("");
   MV /= meanValues.size();
   double DV = 0;
   for (int i = 0; i < (int)meanValues.size(); i++) {
      double meanAbs = abs(meanValues[i]);
      DV += (meanAbs-MV)*(meanAbs-MV);
   }
   DV = sqrt(DV);
   printf("Mean = %.12lf\n", MV);
   printf("Standard deviation = %.12lf\n", DV);
   printf("Standard deviation / Mean = %.12lf\n", DV/MV);
   printf("Mean / Standard deviation = %.12lf\n", MV/DV);
   fftrec(meanValues, false);
   printf("DFT: ");
   for (int i = 0; i < (int)meanValues.size(); i++)
      printf("%.12lf\t", abs(meanValues[i]));
   puts("");
   fflush(stdout);
   */
   QImage spectrumImage = writeSpectrumImage(meanValues);
   spectrumImage.save("!Motion-SPECTRUM.png", "png");
   QLabel* imageLabel = new QLabel();
   imageLabel->setPixmap(QPixmap::fromImage(spectrumImage));
   imageLabel->setFixedSize(spectrumImage.width(), spectrumImage.height());
   imageLabel->show();
}

QImage LaserEnv::writeSpectrumImage(const vector<cc>& row) {
   int len = row.size();
   double mini = 1e10, maxi = -1e10;
   for (int i = 0; i < len; i++) {
      double val = abs(row[i]);
      mini = min(mini, val);
      maxi = max(maxi, val);
   }
   int step = 20;
   int height = 100;
   int width = (len+1)*step;
   QImage retImage = QImage(width, height, QImage::Format_Mono);
   for (int i = 0; i < width; i++)
      for (int j = 0; j < height; j++)
         retImage.setPixel(i, j, 0);
   for (int i = 0; i < width; i++)
      retImage.setPixel(i, height-10, 1);
   for (int i = 0; i < len; i++) {
      double val = (abs(row[i]) - mini) / (maxi - mini);
      int h = val * (height-10);
      for (int j = 0; j < h; j++)
         retImage.setPixel(step*(i+1), height-j-10, 1);
   }
   return retImage;
}

LaserEnv::LaserEnv(QWidget* obj)
   : QWidget(obj), n(256)
{
   circleNumber = 0;
   sz = 10;
   w = 40;
   ray1X = 95;
   ray1Y = 95;
   ray2X = 165;
   ray2Y = 165;
   deviation = 1;
   averageIterations = 5;
   motionSpeed = 4;
   frameX = 150;
   frameY = 150;
   frameHeight = 50;
   frameWidth = 50;
   pointY = 200;
   createGUIControls();
   setupGUIControls();
   //img1.save("initialDFTImage.bmp", "bmp");
   //img2.save("initialRay.bmp", "bmp");
   //img3.save("initialLattice.bmp", "bmp");
}
