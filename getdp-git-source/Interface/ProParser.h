// GetDP - Copyright (C) 1997-2020 P. Dular and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/getdp/getdp/issues.

#ifndef PRO_PARSER_H
#define PRO_PARSER_H

#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <vector>
#include <string>
#include "ListUtils.h"
#include "Message.h"

struct Constant {
  char *Name;
  int Type;
  union {
    double Float;
    char *Char;
    List_T *List;
  } Value;
};

#define VAR_FLOAT         1
#define VAR_LISTOFFLOAT   2
#define VAR_CHAR          3
#define VAR_LISTOFCHAR    4

struct TwoChar { char *char1, *char2; };

// classes for Struct
class Struct {
public:
  Struct() {}
  Struct(int tag,
         std::map<std::string, std::vector<double> > & fopt,
         std::map<std::string, std::vector<std::string> > & copt,
         int member_ValMax) :
    _tag(tag), _member_ValMax(member_ValMax), _fopt(fopt), _copt(copt) {}
  ~Struct() {}

  int append(int tag,
             std::map<std::string, std::vector<double> > & fopt,
             std::map<std::string, std::vector<std::string> > & copt,
             int member_ValMax)
  {
    if (tag >= 0) _tag = tag;
    _member_ValMax = member_ValMax;
    _fopt.insert(fopt.begin(), fopt.end());
    _copt.insert(copt.begin(), copt.end());
    return _tag;
  }

  inline int getTag() const { return _tag; }

  int getMember (std::string & key_member, double & out, int index) const
  {
    std::map<std::string, std::vector<double> >::const_iterator
      it = _fopt.find(key_member);
    if (it != _fopt.end()) {
      if (index < (int)it->second.size()) {
        out = it->second[index]; return 0;
      }
      else {
        out = 0.; return 2; // Error: Index out of range
      }
    }
    else {
      out = 0.; return 1; // Error: Unknown member of Struct
    }
  }

  int getMember (std::string & key_member, const std::string * & out, int index) const
  {
    std::map<std::string, std::vector<std::string> >::const_iterator
      it = _copt.find(key_member);
    if (it != _copt.end()) {
      if (index < (int)it->second.size()) {
        out = &it->second[index]; return 0;
      }
      else {
        out = NULL; return 2; // Error: Index out of range
      }
    }
    else {
      out = NULL; return 1; // Error: Unknown member of Struct
    }
  }

  int getMember_Dim (std::string & key_member, int & out) const
  {
    std::map<std::string, std::vector<double> >::const_iterator
      it = _fopt.find(key_member);
    if (it != _fopt.end()) {
      out = it->second.size(); return 0;
    }
    else {
      std::map<std::string, std::vector<std::string> >::const_iterator
        it = _copt.find(key_member);
      if (it != _copt.end()) {
        out = it->second.size(); return 0;
      }
      else {
        out = 0; return 1; // Error: Unknown member of Struct
      }
    }
  }

  int getMember_String_Dim (std::string & key_member, int & out) const
  {
    std::map<std::string, std::vector<std::string> >::const_iterator
      it = _copt.find(key_member);
    if (it != _copt.end()) {
      out = it->second.size(); return 0;
    }
    else {
      out = 0; return 1; // Error: Unknown member of Struct
    }
  }

  int getMember_Vector (std::string & key_member, const std::vector<double> * & out_vector) const
  {
    std::map<std::string, std::vector<double> >::const_iterator
      it = _fopt.find(key_member);
    if (it != _fopt.end()) {
      out_vector = &it->second; return 0;
    }
    else {
      out_vector = NULL; return 1; // Error: Unknown member of Struct
    }
  }

  int getMember_Vector (std::string & key_member, const std::vector<std::string> * & out_vector) const
  {
    std::map<std::string, std::vector<std::string> >::const_iterator
      it = _copt.find(key_member);
    if (it != _copt.end()) {
      out_vector = &it->second; return 0;
    }
    else {
      out_vector = NULL; return 1; // Error: Unknown member of Struct
    }
  }

  int getMember_ValMax () const { return _member_ValMax; }

  void sprint(std::string & str,
              const std::string & struct_name, const std::string & struct_namespace)
    const
  {
    str = "Struct ";
    if (struct_namespace.size()) str += struct_namespace + "::";
    str += struct_name + " [ ";
    bool flag_comma = false;
    for (std::map<std::string, std::vector<double> >::const_iterator
           it_attrib = _fopt.begin();
         it_attrib != _fopt.end(); ++it_attrib ) {
      if (!flag_comma && it_attrib != _fopt.begin()) flag_comma = true;
      if (flag_comma) str += ", ";
      str += it_attrib->first + " ";
      char tmp[32];
      if (it_attrib->second.size() > 1) str += "{ ";
      for (unsigned int i = 0; i < it_attrib->second.size(); i++) {
        if (i) str += ", ";
        sprintf(tmp, "%g", it_attrib->second[i]); str += tmp;
      }
      if (it_attrib->second.size() > 1) str += "}";
    }
    for (std::map<std::string, std::vector<std::string> >::const_iterator
           it_attrib = _copt.begin();
         it_attrib != _copt.end(); ++it_attrib ) {
      if (!flag_comma && it_attrib != _copt.begin()) flag_comma = true;
      if (flag_comma) str += ", ";
      str += it_attrib->first + " ";
      if (it_attrib->second.size() > 1) str += "Str[{ ";
      for (unsigned int i = 0; i < it_attrib->second.size(); i++) {
        if (i) str += ", ";
        str += "\"" + it_attrib->second[i] + "\"";
      }
      if (it_attrib->second.size() > 1) str += "}]";

    }
    str += " ];\n";
  }

private:
  int _tag, _member_ValMax;
  std::map<std::string, std::vector<double> > _fopt;
  std::map<std::string, std::vector<std::string> > _copt;
};


template <class K, class T>
class Map {
public:
  Map() {}
  ~Map() {}

  T * Find(K key)
  {
    typename std::map<K, T>::iterator it;
    if ( (it = _map.find(key)) != _map.end() ) return &it->second;
    else return NULL;
  }

  const T * Find(K key) const
  {
    typename std::map<K, T>::const_iterator it;
    if ( (it = _map.find(key)) != _map.end() ) return &it->second;
    else return NULL;
  }

  inline T & operator[] (K key) { return _map[key]; }
  inline std::map<K, T> & get() { return _map; }
  inline const std::map<K, T> & get() const { return _map; }
  inline int count (const std::string key) const { return _map.count(key); }
  inline int size () const { return _map.size(); }
  void clear() { _map.clear(); }

public:
  std::map<K, T> _map;
};


typedef std::map<std::string, Struct> Map_string_Struct;

class Structs : public Map<std::string, Struct> {
public:
  Structs() { _max_tag = 0; }
  ~Structs() {}

  int defStruct(std::string & struct_name,
                std::map<std::string, std::vector<double> > & fopt,
                std::map<std::string, std::vector<std::string> > & copt,
                int member_ValMax, bool append = false)
  {
    int tag;
    std::map<std::string, std::vector<double> >::const_iterator it = fopt.find("Tag");
    if (it != fopt.end()) {
      tag = (int)it->second[0]; // Tag forced
      _max_tag = std::max(_max_tag, tag);
    }
    else {
      tag = (!append)? ++_max_tag : -1; // Tag auto
      if (!append) fopt["Tag"].push_back((double)tag);
    }
    if (!append)
      (*this)[struct_name] = Struct(tag, fopt, copt, member_ValMax);
    else
      (*this)[struct_name].append(tag, fopt, copt, member_ValMax);
    return tag;
  }

  int get_key_struct_from_tag(int tag, const std::string * & key_struct) const
  {
    Map_string_Struct::const_iterator it_st;
    for (it_st = this->get().begin(); it_st != this->get().end(); ++it_st )
      if (it_st->second.getTag() == tag) break;
    if (it_st == this->get().end()) return 2; // 2: Error: Unknown Struct
    key_struct = &it_st->first;
    return 0; // 0: no error
  }

  void sprint(std::vector<std::string> & strs_out, const std::string & struct_namespace) const
  {
    std::string str;
    for (Map_string_Struct::const_iterator it_st = this->get().begin();
         it_st != this->get().end(); ++it_st ) {
      it_st->second.sprint(str, it_st->first, struct_namespace);
      strs_out.insert(strs_out.end(), str);
    }
  }

private:
  int _max_tag;
};


typedef std::map<std::string, Structs> Map_string_Structs;

class NameSpaces : public Map<std::string, Structs> {
public:
  NameSpaces() {}
  ~NameSpaces() {}

  int defStruct(std::string & key_namespace, std::string & key_name,
                std::map<std::string, std::vector<double> > & fopt,
                std::map<std::string, std::vector<std::string> > & copt,
                int & tag_out, int member_ValMax, bool append = false)
  {
    Structs * structs_P = &(*this)[key_namespace];
    if (structs_P->count(key_name)) {
      if (!append) {
        tag_out = (*structs_P)[key_name].getTag();
        return 1; // 1: Error: Redefinition of Struct
      }
    }
    else if (append) append = false; // non-existing Struct
    tag_out = structs_P->defStruct(key_name, fopt, copt, member_ValMax, append);
    return 0; // 0: no error
  }

  int getTag(std::string & key_namespace, std::string & key_name,
             double & out) const
  {
    const Structs * structs_P = this->Find(key_namespace);
    const Struct * struct_P = (structs_P)? structs_P->Find(key_name) : NULL;
    if (structs_P && struct_P) {
      out = (double)struct_P->getTag();
    }
    else  {
      out = 0.; return 1; // 1: Error: Unknown Struct
    }
    return 0; // 0: no error
  }

  int getMember(std::string & key_namespace, std::string & key_name,
                std::string & key_member, double & out, int index = 0) const {

    const Structs * structs_P = this->Find(key_namespace);
    const Struct * struct_P = (structs_P)? structs_P->Find(key_name) : NULL;
    if (structs_P && struct_P) {
      switch (struct_P->getMember(key_member, out, index)) {
      case 0:
        break;
      case 1:
        out = 0.; return 2; // 2: Error: Unknown member of Struct
        break;
      case 2:
        out = 0.; return 3; // 3: // Error: Index out of range
        break;
      }
    }
    else  {
      out = 0.; return 1; // 1: Error: Unknown Struct
    }
    return 0; // 0: no error
  }

  int getMember(std::string & key_namespace, std::string & key_name,
                std::string & key_member, const std::string * & out, int index = 0) const
  {

    const Structs * structs_P = this->Find(key_namespace);
    const Struct * struct_P = (structs_P)? structs_P->Find(key_name) : NULL;
    if (structs_P && struct_P) {
      switch (struct_P->getMember(key_member, out, index)) {
      case 0:
        break;
      case 1:
        out = NULL; return 2; // 2: Error: Unknown member of Struct
        break;
      case 2:
        out = NULL; return 3; // 3: // Error: Index out of range
        break;
      }
    }
    else  {
      out = NULL; return 1; // 1: Error: Unknown Struct
    }
    return 0; // 0: no error
  }

  int getMember_Dim(std::string & key_namespace, std::string & key_name,
                    std::string & key_member, int & out) const {

    const Structs * structs_P = this->Find(key_namespace);
    const Struct * struct_P = (structs_P)? structs_P->Find(key_name) : NULL;
    if (structs_P && struct_P) {
      switch (struct_P->getMember_Dim(key_member, out)) {
      case 0:
        break;
      case 1:
        out = 0; return 2; // 2: Error: Unknown member of Struct
        break;
      }
    }
    else  {
      out = 0; return 1; // 1: Error: Unknown Struct
    }
    return 0; // 0: no error
  }

  int getMember_Vector(std::string & key_namespace, std::string & key_name,
                       std::string & key_member, const std::vector<double> * & out_vector) const {

    const Structs * structs_P = this->Find(key_namespace);
    const Struct * struct_P = (structs_P)? structs_P->Find(key_name) : NULL;
    if (structs_P && struct_P) {
      switch (struct_P->getMember_Vector(key_member, out_vector)) {
      case 0:
        break;
      case 1:
        out_vector = NULL; return 2; // 2: Error: Unknown member of Struct
        break;
      }
    }
    else  {
      out_vector = NULL; return 1; // 1: Error: Unknown Struct
    }
    return 0; // 0: no error
  }

  int getMember_Vector(std::string & key_namespace, std::string & key_name,
                       std::string & key_member, const std::vector<std::string> * & out_vector) const {

    const Structs * structs_P = this->Find(key_namespace);
    const Struct * struct_P = (structs_P)? structs_P->Find(key_name) : NULL;
    if (structs_P && struct_P) {
      switch (struct_P->getMember_Vector(key_member, out_vector)) {
      case 0:
        break;
      case 1:
        out_vector = NULL; return 2; // 2: Error: Unknown member of Struct
        break;
      }
    }
    else  {
      out_vector = NULL; return 1; // 1: Error: Unknown Struct
    }
    return 0; // 0: no error
  }

  int get_key_struct_from_tag(std::string & key_namespace,
                              int tag, const std::string * & key_struct) const
  {
    const Structs * structs_P = this->Find(key_namespace);
    if (structs_P != NULL)
      return structs_P->get_key_struct_from_tag(tag, key_struct);
    else return 1; // 1: Error: Unknown NameSpace
  }

  int getMember_ValMax(std::string & key_namespace, std::string & key_name)
  {
    const Structs * structs_P = this->Find(key_namespace);
    const Struct * struct_P = (structs_P)? structs_P->Find(key_name) : NULL;
    return (structs_P && struct_P)? struct_P->getMember_ValMax() : -1;
  }

  void sprint(std::vector<std::string> & strs_out) const
  {
    std::vector<std::string> strs;
    for (Map_string_Structs::const_iterator it_ns = this->get().begin();
         it_ns != this->get().end(); ++it_ns ) {
      strs.clear();
      it_ns->second.sprint(strs, it_ns->first);
      strs_out.insert(strs_out.end(), strs.begin(), strs.end());
    }
  }
};


extern FILE *getdp_yyin;
extern std::string getdp_yyname;
extern char getdp_yyincludename[256];
extern long int getdp_yylinenum;
extern int getdp_yycolnum;
extern int getdp_yyincludenum;
extern int level_include;
extern int getdp_yyerrorlevel;
extern std::map<std::string, std::vector<double> > CommandLineNumbers;
extern std::map<std::string, std::vector<std::string> > CommandLineStrings;
extern std::map<std::string, std::vector<double> > GetDPNumbers;
extern std::map<std::string, std::vector<std::string> > GetDPStrings;
extern std::map<std::string, std::map<int, std::vector<double> > > GetDPNumbersMap;

int getdp_yyparse();
void getdp_yyrestart(FILE*);
void Free_ParserVariables();
char *strSave(const char *string);
char *strEmpty();
void cStyleComments();
void cxxStyleComments();
void parseString(char endchar);
void skipUntil(const char *skip, const char *until);
void skipUntil_test(const char *skip, const char *until,
                    const char *until2, int l_until2_sub, int *type_until2);
void Print_Constants();
void Print_Struct();
int  Print_ListOfDouble(char *format, List_T *list, char *buffer);
Constant *Get_ParserConstant(char *name);

#endif
