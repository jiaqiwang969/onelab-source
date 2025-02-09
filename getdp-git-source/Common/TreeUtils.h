// GetDP - Copyright (C) 1997-2020 P. Dular and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/getdp/getdp/issues.

#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include "avl.h"
#include "ListUtils.h"

typedef struct {
  int size;
  avl_tree *root;
} Tree_T;

Tree_T *Tree_Create(int size, int (*fcmp)(const void *a, const void *b));
void    Tree_Delete(Tree_T *Tree);
void   *Tree_Add(Tree_T *tree, void *data);
int     Tree_Nbr(Tree_T *Tree);
int     Tree_Insert(Tree_T *Tree, void *data);
int     Tree_Replace(Tree_T *Tree, void *data);
int     Tree_Search(Tree_T *Tree, void *data);
int     Tree_Query(Tree_T *Tree, void *data);
void   *Tree_PQuery(Tree_T *Tree, void *data);
int     Tree_Suppress(Tree_T *Tree, void *data);
int     Tree_Size(Tree_T *tree);
void    Tree_Action(Tree_T *tree, void (*action) (void *data, void *dummy));
List_T *Tree2List(Tree_T *pTree);

#endif
