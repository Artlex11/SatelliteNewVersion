#include "Links.h"

// ���������� ������ ���������� ����� (�������� �� ������)
//void Links::addLink(LinkData& link) {
//    links.push_back(link);
//}

// ���������� ������ ���������� ����� � �������������� �����������
void Links::addLink(LinkData&& link) {
    links_vec.push_back(std::move(link));
}

// ���������� ������ ��������� ���� ������
std::vector<LinkData>& Links::getLinks() {
    return links_vec;
}