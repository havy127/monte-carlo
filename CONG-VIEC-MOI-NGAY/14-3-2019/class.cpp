/**
 * User: HaVy HUYNH
 * Date: 03/14/2019
 * Time: 9:18 AM
 */
#include <iostream>
#include <string>
using namespace std;

int BienX;

//Level 1
class LopHocLapTrinh{
  public:
    void lopHocLapTrinhCPlus() {
      std::cout << "Chao cac ban !!!!" << '\n';
    }
};
// Level 2
class User{
  private:
    string userName;
  public:
    void setName(string Name) {
      userName = Name;
    }
    string getName(){
      return "Usernam la: " + userName + '\n';
    }
};


int main(int argc, char const *argv[]) {
  //level 1
  // LopHocLapTrinh LopHocLapTrinhObj;
  // LopHocLapTrinhObj.lopHocLapTrinhCPlus();
  //level 2
  User username;
  username.setName("HUYNH");
  std::cout << "Chao ban !!!!" <<username.getName() <<'\n';
  return 0;
}
