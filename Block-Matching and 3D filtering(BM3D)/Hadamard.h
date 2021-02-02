#include <iostream>
#include <vector>

//生成哈达玛矩阵的函数
int CreateHadmard(int i, int j)
{
    long k, temp, result = 0;
    temp = i & j;
    for (k = 0; k < 32; k++)
    {
        result = result + (temp >> k) & 1;
    }
    if (result % 2 == 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

void Hadmard(std::vector<int> &Hadamard,int M,int N)
{

    int i, j;

    printf("生成哈达玛矩阵\n");
    for (i = 0; i < M*N; i++)
    {   
        //std::vector<int> temp;
        //Hadamard.push_back(temp);
		/*for (j = 0; j < N; j++)
		{*/
            Hadamard[i].push_back(CreateHadmard(i, j));
        //}
    }
}