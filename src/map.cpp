#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
struct pairhash{
    size_t operator()(const std::pair<int, int> &p) const {
        return
            std::hash<unsigned int>()(p.first) ^
            std::hash<unsigned int>()(p.second);
    }
};
typedef std::unordered_map<std::pair<int, int>, std::vector<int>*, pairhash> hbond_map;
typedef std::unordered_map<int, std::vector<int>*> hbond_int_map;
typedef std::unordered_map<std::pair<int, int>, std::vector<int>*, pairhash>::iterator hbond_map_itr;
typedef std::unordered_map<int, std::vector<int>*>::iterator hbond_int_map_itr;

int main ()
{
    using namespace std;
    unordered_map<int,double> mymap = {
        {5,5.4},
        {4,6.1},
        {2,5.9} };
    hbond_int_map test;
    hbond_map* hbonds;
    pair<int, int> hbond_id(540,30);
    int* newvector= vector<int>*();
    hbonds->insert(make_pair(hbond_id,vector<int>*[(0)]));
    printf("hbond_int_map test %d\n",test.count(0));

    int input=6;

    unordered_map<int,double>::iterator got = mymap.find (input);
    cout << mymap.count(input) << endl;

    //got->second=6.0;

    //cout << got->first << " is " << got->second;

    //cout << endl;

    return 0;
}
