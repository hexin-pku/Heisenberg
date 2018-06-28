#ifndef basic_Tools_H
#define basic_Tools_H

#include <string>

using namespace std;

template <class T>
int getArrayLen(T& array)
{
return (sizeof(array) / sizeof(array[0]));
};

string& trim(string &s)   
{  
    if (s.empty())   
    {  
        return s;  
    }  
    s.erase(0,s.find_first_not_of(" "));
    s.erase(s.find_last_not_of("\r") + 1);  
    s.erase(s.find_last_not_of(" ") + 1);
    return s;  
}

string& replace_recursive(string& str, const string& old_value, const string& new_value)     
{
    while(true)   {
        string::size_type pos(0);     
        if(   (pos=str.find(old_value))!=string::npos   )     
            str.replace(pos,old_value.length(),new_value);     
        else break;     
    }     
    return str;     
}

string& replace_distinct(string& str, const string& old_value, const string& new_value)     
{
    for(string::size_type pos(0); pos!=string::npos; pos+=new_value.length())
    {     
        if( (pos=str.find(old_value,pos))!=string::npos )     
            str.replace(pos,old_value.length(),new_value);     
        else break;     
    }     
    return str;     
}

#endif


