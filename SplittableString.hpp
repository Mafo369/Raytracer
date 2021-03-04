#ifndef SPLITTABLE_STRING
#define SPLITTABLE_STRING


#include <string>
#include <vector>
#include <iostream>

class SplittableString : public std::string
{
public : 
    SplittableString() : std::string(){}
    SplittableString(const std::string& str) : std::string(str){}
    SplittableString(const char* str) : std::string(str){}

    std::vector<SplittableString> split(const std::string& sequence) const
    {
        std::vector<SplittableString> result;

        size_t i = 0;
        size_t split_position;
        while (i < size() && (split_position = find(sequence, i)) != std::string::npos)
        {
            result.push_back(substr(i, split_position - i));
            i = split_position + sequence.length();
        }
        
        if (i < size())
        {
            result.push_back(substr(i));
        }

        return result;
    }

};



// int main(void)
// {
//     SplittableString str = "f 102/347/121 118/348/121 103/349/121";
//     for (const SplittableString& s: str.split(" "))
//     {
//         for (const SplittableString& n : s.split("/"))
//         {
//             std::cout << n << "&";
//         }

//         std::cout << "\n" << std::endl;
//     }
//     return 0;
// }



#endif