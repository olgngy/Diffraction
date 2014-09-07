#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPaintEvent>
#include <QLineEdit>
#include <QPushButton>

class MainWindow : public QMainWindow {
    Q_OBJECT

    int WINDOW_WIDTH;
    int WINDOW_HEIGHT;
    int WINDOW_IDENT;
    double CIRCLE_CENTER_X;
    double CIRCLE_CENTER_Y;
    double CIRCLE_RADIUS;
    double RAY_Y;
    double N1;
    double N2;

    double L;
    double Phi;

    QLineEdit* editCircleCenterX;
    QLineEdit* editCircleCenterY;
    QLineEdit* editCircleRadius;
    QLineEdit* editRayY;
    QLineEdit* editN1;
    QLineEdit* editN2;
    QPushButton* updateButton;

    void installUserInterface();

    void paintEvent(QPaintEvent* paintEvent);
    void drawPlot(QPainter* painter);
    void drawRay(QPainter* painter);
    void drawUserInterface(QPainter* painter);

private slots:
    void updateButtonClicked();
    
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();  
};

#endif // MAINWINDOW_H
