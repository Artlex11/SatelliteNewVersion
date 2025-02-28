#include "Links.h"

// Реализация метода добавления связи (передача по ссылке)
//void Links::addLink(LinkData& link) {
//    links.push_back(link);
//}

// Реализация метода добавления связи с использованием перемещения
void Links::addLink(LinkData&& link) {
    links_vec.push_back(std::move(link));
}

// Реализация метода получения всех связей
std::vector<LinkData>& Links::getLinks() {
    return links_vec;
}