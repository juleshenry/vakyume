import os
import shutil

def test_case_sensitivity():
    test_dir = "test_case_sensitivity"
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir)
    
    file_P = os.path.join(test_dir, "eqn_1_13b__P.py")
    file_p_a = os.path.join(test_dir, "eqn_1_13b__p_a.py")
    
    with open(file_P, "w") as f:
        f.write("P = 1")
    
    print(f"Created {file_P}")
    
    with open(file_p_a, "w") as f:
        f.write("p_a = 2")
    
    print(f"Created {file_p_a}")
    
    files = sorted(os.listdir(test_dir))
    print(f"Files in directory: {files}")

if __name__ == "__main__":
    test_case_sensitivity()
