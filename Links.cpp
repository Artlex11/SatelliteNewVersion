#include "links.h"


// Реализация метода добавления связи с использованием перемещения
void Links::addLink(LinkData&& link) {
    links.emplace_back(std::move(link));
}

// Реализация метода получения всех связей
std::vector<LinkData>& Links::getLinks() {
    return links;
}