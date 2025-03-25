#include "links.h"


// ���������� ������ ���������� ����� � �������������� �����������
void Links::addLink(LinkData&& link) {
    links.emplace_back(std::move(link));
}

// ���������� ������ ��������� ���� ������
std::vector<LinkData>& Links::getLinks() {
    return links;
}