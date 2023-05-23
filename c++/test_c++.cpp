#include <iostream>
#include <cassert>

// Function to test
int add(int a, int b) {
  return a + b;
}

// Test cases
void test_add() {
  assert(add(1, 2) == 3);
  assert(add(0, 0) == 0);
  assert(add(-1, 1) == 0);
  std::cout << "All test cases for add() passed." << std::endl;
}

// Test runner
int main() {
  test_add();
  return 0;
}