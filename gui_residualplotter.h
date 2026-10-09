#pragma once
#include <QWidget>
#include <vector>

class ResidualPlotter : public QWidget
{
    Q_OBJECT
public:
    explicit ResidualPlotter(QWidget* parent = nullptr);
    void setData(const std::vector<double>& residuals);

protected:
    void paintEvent(QPaintEvent*) override;

private:
    std::vector<double> m_data;
};
