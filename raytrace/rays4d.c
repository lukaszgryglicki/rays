#define DIM_4D
#define REAL __REAL_64bit
#include <string.h>
#include <math.h>
#include "rayslib.h"

REAL slimit_min, slimit_max, slimit_step;
int vx_ndiv, vx_whendiv;
BTree* btree;
BList* boxes;
Triangle* g_ts;
int loaded;
Vertex light, observer;
Screen screen;
IdxTree* itree;
REAL lookz;
int line_idx;
int (*intersection)(Triangle*, Ray*, Vertex*, int*, IdxTree*);
void (*get_normal)(Vector*, Triangle*, Vertex*, Material*, Material*, Material*, TexCoord*);
int (*intersection_itree)(Ray* r, IdxTree* i_tree);

void get_normal_4d(Vector*, Triangle*, Vertex*, Material*, Material*, Material*, TexCoord*);
int intersection_itree_4d(Ray* r, IdxTree* i_tree);
int intersection_4d(Triangle*, Ray*, Vertex*, int*, IdxTree*);

void load_scene(Triangle** ts, char* scenefile)
{
 FILE* f;
 char texDir[1024];
 char str[1024];
 char str2[1024];
 REAL dd;
 /*int texid;*/
 int nr,sx,sy,n,i,i1,i2,ii1,ii2;
 long pos;
 int gtrans;
 REAL **lm, **lmn;
 int *t_use;
 char** tmap;
 ListTransform* ltemp;
 MaterialTransform omt;
 args_fprintf(stdout, "Reading scene definition...\n");
 ltemp = NULL;
 tlist = NULL;
 tmap = NULL;
 trans_used = 0;
 vlight = 0;
 lookz = .333333;
 world = matrix(4);
 I_matrix(world, 4);
 worldn = matrix(4);
 I_matrix(worldn, 4);
 f = fopen(scenefile,"rb");	/* is ok with r? */
 if (!f) spanic(HERE, "load_scene: cant open file '%s'", scenefile);
 if (is_binary(f)) 
 { 
     load_binary_scene(f, ts); 
     return; 
 }
 nr = v_fscanf(f, "Screen: (%d,%d)\n", &sx, &sy);
 if (nr != 2) spanic(HERE, "load_scene: cant read screen def", HERE);
 if (sx <= 0 || sy <= 0) spanic(HERE, "load_scene: bad screen def: %dx%d", sx, sy);
 if (ovr_x > 0) sx = ovr_x;
 if (ovr_y > 0) sy = ovr_y;
 if (double_res) { sx *= 2; sy *= 2; }
 pos = ftell(f);
 nr = v_fscanf(f, "%s", str);
 if (nr != 1) spanic(HERE, "load_scene: cant read background option or other", HERE);
 if (!strcmp(str, "Background:"))
   {
    nr = v_fscanf(f, " %d\n", &bkgnd);
    if (nr != 1) spanic(HERE, "load_scene: cant read background definition", HERE);
    if (bkgnd < 0.) spanic(HERE, "bkgnd lower than 0", HERE);
   }
 else 
   {
    fseek(f, pos, SEEK_SET);
    bkgnd = 0;
   }
 pos = ftell(f);
 nr = v_fscanf(f, "%s", str);
 if (nr != 1) spanic(HERE, "load_scene: cant read global norm or other option", HERE);
 if (!strcmp(str, "NormalDistorber:"))
   {
    nr = v_fscanf(f, " %~%%\n", &dd);
    if (nr != 1) spanic(HERE, "load_scene: cant read normal distorber definition", HERE);
    if (dd < 0.) dd = 0.;
    if (!ovr_no) global_dist = dd/100.;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f, "%s", str);
 if (nr != 1) spanic(HERE, "load_scene: cant read shadow or recurse defs", HERE);
 if (!strcmp(str, "MinShadow:"))
   {
    nr = v_fscanf(f, " %~\n", &dd);
    if (nr != 1) spanic(HERE, "load_scene: cant read minshadow", HERE);
    if (dd < 0.) dd = 0.;
    if (dd > 1.) dd = 1.;
    if (!ovr_m) minshadow = 1. - dd;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f, "%s", str);
 if (nr != 1) spanic(HERE, "load_scene: cant read shadow or recurse defs", HERE);
 if (!strcmp(str, "MaxShadow:"))
   {
    nr = v_fscanf(f, " %~\n", &dd);
    if (nr != 1) spanic(HERE, "load_scene: cant read maxshadow", HERE);
    if (dd < 0.) dd = 0.;
    if (dd > 1.) dd = 1.;
    if (!ovr_ms) maxshadow = 1. - dd;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f, "%s", str);
 if (nr != 1) spanic(HERE, "load_scene: cant read shadow or recurse defs", HERE);
 if (!strcmp(str, "Ambient:"))
   {
    nr = v_fscanf(f, " %~\n", &dd);
    if (nr != 1) spanic(HERE, "load_scene: cant read ambient", HERE);
    if (dd < 0.) dd = 0.;
    if (dd > 1.) dd = 1.;
    if (!ovr_a) ambient = dd;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f, "%s", str);
 if (nr != 1) spanic(HERE, "load_scene: cant read shadow or recurse defs", HERE);
 if (!strcmp(str, "MaxRecurse:"))
   {
    nr = v_fscanf(f, " %d\n", &i);
    if (nr != 1) spanic(HERE, "load_scene: cant read maxrecurse", HERE);
    if (i < 0) i = 0;
    if (!ovr_r) max_rec = i;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f, "%s", str);
 if (nr != 1) spanic(HERE, "load_scene: cant read backup or observer def", HERE);
 if (!strcmp(str, "Backup:"))
   {
    nr = v_fscanf(f, " %d\n", &i);
    if (nr != 1) spanic(HERE, "load_scene: cant read backup definition", HERE);
    if (i < 8) i = 8;
    if (!ovr_b) bkup = i;
   }
 else fseek(f, pos, SEEK_SET);
 init_screen(&screen, sx, sy);
 v_fscanf(f, "Observer: ");
 if (!read_vertex(f, &observer)) spanic(HERE, "load_scene: cant read observer def", HERE);
 v_fscanf(f,"\n");
 pos = ftell(f);
 nr = v_fscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "ObserverTransform:"))
       {
/*	   args_fprintf(stdout, "Observer transform.\n");*/
	lm = matrix(4);
	lmn = matrix(4);
	I_matrix(lm, 4);
	I_matrix(lmn, 4);
	clear_material_transform(&omt);
 	read_transformation_long(f, &lm, &lmn, &i, &omt);
	transform_observer(&observer, lm, lmn);
	free_matrix(lm, 4);
	free_matrix(lmn, 4);
	lm = lmn = NULL;
       }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "LookZ:"))
   {
    nr = v_fscanf(f, "%~%%\n", &lookz);
    if (nr != 1) spanic(HERE, "load_scene: cant read lookZ def", HERE);
    lookz /= 100.;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "Light:"))
   {
    if (!read_vertex(f, &light)) spanic(HERE, "load_scene: cant read light vertex def", HERE);
   }
 else if (nr == 1 && !strcmp(str, "VLight:"))
   {
    if (!read_vector(f, &light)) spanic(HERE, "load_scene: cant read light vector def", HERE);
    vlight = 1;
    if (nearly_equal(length(&light), 0., precision)) spanic(HERE, "load_scene: vector light: length too close 0\n", HERE);
    normalize(&light);
   }
 else spanic(HERE, "load_scene: cant read light def", HERE);
 pos = ftell(f);
 nr = v_fscanf(f,"%s %s\n", str, str2);
 if (nr == 2 && !strcmp(str, "LightColor:"))
       {
        if (!l_global)
		  {
		   set_light_col(str2);
		   check_light();
		  }
       }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "LightTransform:"))
       {
	lm = matrix(4);
	lmn = matrix(4);
	I_matrix(lm, 4);
	I_matrix(lmn, 4);
	clear_material_transform(&omt);
 	read_transformation_long(f, &lm, &lmn, &i, &omt);
 	/*print_matrix(lm, 4);
	print_matrix(lmn, 4);
	args_fprintf(stdout, "done.\n");*/
	transform_light(&light, lm, lmn, vlight);
/* args_fprintf(stdout, "(%~,%~,%~)\n", light.x, light.y, light.z);*/
	free_matrix(lm, 4);
	free_matrix(lmn, 4);
	lm = lmn = NULL;
       }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = v_fscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "WorldTransform:"))
       {
		clear_material_transform(&omt);
		read_transformation_long(f, &world, &worldn, &trans_used, &omt);
        if (dd > 0. && !ovr_no && global_dist <= 0.)
	  {
	   global_dist = dd;
/*	   args_fprintf(stdout, "world_transform -> set global_dist to: %~\n", global_dist);*/
	  }
       }
 else fseek(f, pos, SEEK_SET);
 if (awtmat) 
 {
     awt_apply(&trans_used);
     print_matrix(world, 4);
 }
 nr = v_fscanf(f, "TexDirectory: %s\n", texDir);
 if (nr != 1) spanic(HERE, "load_scene: cant read texture directory", HERE);
 nr = v_fscanf(f, "NumTextures: %d\n", &nTex);
 if (nr != 1) spanic(HERE, "load_scene: cant read num textures", HERE);
 if (nTex < 0) spanic(HERE, "load_scene: negative textures count", HERE);
 if (bkgnd < 0 || bkgnd > nTex) spanic(HERE, "bkgnd index out of range", HERE);
 create_textures(texDir);
 if (texture)
  {
   t_use = (int*)malloc(nTex*sizeof(int));
   for (i=0;i<nTex;i++) t_use[i] = 0;
   if (bkgnd > 0) t_use[bkgnd-1] ++;
  }
 else t_use = NULL;
 pos = ftell(f);
 nr = v_fscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "TextureMapping:"))
   {
    v_fscanf(f,"{\n");
    tmap = (char**)malloc(nTex*sizeof(char*));
    for (i=0;i<nTex;i++) tmap[i] = NULL;
    while (1)
      {
       pos = ftell(f);
       nr = v_fscanf(f, "%s\n", str);
       if (nr == 1 && !strcmp(str,"}")) break;
       else fseek(f, pos, SEEK_SET);
       nr = v_fscanf(f, "%d:", &i);
       if (nr != 1) spanic(HERE, "read_triangle: cant read mapping idx", HERE);
       if (i < 1 || i > nTex) spanic(HERE, "read_triangle: bad mapping idx", HERE);
       i--;
       if (tmap[i]) spanic(HERE, "read_trianmgle: mapping index already used", HERE);
       tmap[i] = (char*)malloc(512*sizeof(char));
       nr = v_fscanf(f,"%s\n", tmap[i]);
       if (nr != 1) spanic(HERE, "load_scene: cant read texture mapping name\n", HERE);
      }
   }
 else fseek(f, pos, SEEK_SET);
 nr = v_fscanf(f, "nTriangles: %d\n", &n);
 if (nr != 1) spanic(HERE, "load_scene: cant read num triangles", HERE);
 if (n < 0) spanic(HERE, "load_scene: negative triangles count", HERE);
 nNURBS = 0;
 pos = ftell(f);
 nr = v_fscanf(f, "nNURBS: %d\n", &i);
 if (nr != 1) fseek(f, pos, SEEK_SET);
 else
   {
    if (i < 0) spanic(HERE, "load_scene: negative NURBS count", HERE);
/*    args_fprintf(stdout, "READING NURBS SURFACE!\n");*/
    nTriNURBS = read_nurbs_from_file(f,i,n);
    args_fprintf(stdout, "Adding %d triangles as pretriangulation of surface..\n", nTriNURBS);
    n += nTriNURBS;
   }
 if (n > 0) *ts = (Triangle*)malloc(n*sizeof(Triangle));
 if (n <= 0 && nNURBS <= 0) spanic(HERE, "nothing to RT", HERE);
 if (nNURBS > 0)
   {
    args_fprintf(stdout, "Starting triangulation, factors: (%1.3~,%1.3~,%1.3~,%1.3~,%1.3~)\n"
	    , tri_pow, tri_pow_edge, tri_pow_edge2, tri_pow_middle, tri_pow_middle2);
    for (i=0;i<nNURBS;i++) triangulate_nurbs(&g_nurbs[i], *ts);
   }
 gtrans = 1;
 while (gtrans)
   {
    pos = ftell(f);
    nr = v_fscanf(f,"%s\n", str);
    if (nr == 1 && !strcmp(str, "ListTransform:"))
       {
	nr = v_fscanf(f, "[%d,%d]", &i1, &i2);
	if (nr != 2) spanic(HERE, "load_scene: bad list transformation (bad indices format)", HERE);
	lm = matrix(4);
	I_matrix(lm, 4);
	lmn = matrix(4);
	I_matrix(lmn, 4);
	
	
	ltemp = (ListTransform*)malloc(sizeof(ListTransform));

	clear_material_transform(&(ltemp->mt));
	read_transformation_long(f, &lm, &lmn, &i, &(ltemp->mt));
	ltemp->M = lm;
	ltemp->MN = lmn;
	ltemp->i1 = i1;
	ltemp->i2 = i2;
	ltemp->ni1 = ltemp->ni2 = -1;
	/*if (dd < 0.) ltemp->n_dist = -1.;
	else ltemp->n_dist = dd;
	if (texid < 0.) ltemp->t_id = -1;
	else ltemp->t_id = texid;*/
	ptrlist_add(&tlist, (void*)ltemp);
       }
    else if (nr == 1 && !strcmp(str, "NURBSTransform:"))
       {
	nr = v_fscanf(f, "[%d,%d]", &i1, &i2);
	if (nr != 2) spanic(HERE, "load_scene: bad nurbs transformation (bad indices format)", HERE);
	ii1 = i1;
	ii2 = i2;
	indices_nurbs2lt(&i1, &i2);
/*	args_fprintf(stdout, "Computed: (%d,%d)\n", i1, i2);*/
	lm = matrix(4);
	I_matrix(lm, 4);
	lmn = matrix(4);
	I_matrix(lmn, 4);
	
	ltemp = (ListTransform*)malloc(sizeof(ListTransform));
	clear_material_transform(&(ltemp->mt));
	read_transformation_long(f, &lm, &lmn, &i, &(ltemp->mt));

	ltemp->M = lm;
	ltemp->MN = lmn;
	ltemp->i1 = i1;
	ltemp->i2 = i2;
	ltemp->ni1 = ii1;
	ltemp->ni2 = ii2;
	/*ltemp->t_id = texid;
	if (dd < 0.) ltemp->n_dist = -1.;
	else ltemp->n_dist = dd;
	if (texid < 0.) ltemp->t_id = -1;
	else ltemp->t_id = texid;*/
	ptrlist_add(&tlist, (void*)ltemp);
       }
   else { fseek(f, pos, SEEK_SET); gtrans = 0; }
   }
 nTriangles = n;
 args_fprintf(stdout, "Reading %d triangles...",nTriangles);
 fflush(stdout);
 for (i=0;i<n-nTriNURBS;i++)
   {
    if (!read_triangle(f, &((*ts)[i]), &i, t_use, tmap, *ts, 1)) spanic(HERE, "load_scene: general triangle read failure: i=%d", i);
/*    args_fprintf(stdout, "(%~,%~,%~)\n", (*ts)[0].a.x, (*ts)[0].a.y, (*ts)[0].a.z);*/
/*    args_fprintf(stdout, "read %d triangles\n", i);*/
   }
 for (i=n-nTriNURBS;i<n;i++)
   {
    if (!read_triangle(f, &((*ts)[i]), &i, t_use, tmap, *ts, 0)) spanic(HERE, "load_scene: general triangle finalize failure, NURBS triangulated, i=%d", i);
   }
 args_fprintf(stdout, "\n");
 update_list_transform(*ts, tlist);
 for (i=0;i<nNURBS;i++)
   {
    list_transform_nurbs(&g_nurbs[i], i);
   }
 if (tmap)
   {
    for (i=0;i<nTex;i++) if (tmap[i]) free(tmap[i]);
    free(tmap);
    tmap = NULL;
   }
 fclose(f);
 check_triangles(*ts);
 if (t_disabled) pernament_disable_texturing(t_use, *ts, nTex, nTriangles);
 postprocess_textures(t_use, *ts, texDir);
 strcpy(textureDir, texDir);
 if (want_bin) 
  {
   save_binary_scene(scenefile, *ts);
/*   save_computed_scene(scenefile, *ts);*/
  }
 if (t_use) free(t_use);
 args_fprintf(stdout, "\nFinally %d triangles.\n", nTriangles);
 debug(HERE, "Finally %d triangles.\n", nTriangles);
 /*for (i=0;i<nTriangles;i++)
 {
  debug(HERE, "%d) %~\n", i, (*ts)[i].a.z);
 }*/
}

void compute_slimits(int n)
{
/* n /= 2;*/
/* args_fprintf(stdout, "nnnnnn = %d\n", n);*/
 if (n < 1000)
   {
    slimit_max = n/2;
    slimit_min = n/2;
    slimit_step = 0.0;
   }
 else if (n >= 1000 && n < 2000)
   {
    slimit_max = n/4;
    slimit_min = n/8;
    slimit_step = (4.*(slimit_max - slimit_min))/n;
   }
 else if (n >= 2000 && n < 4000)
   {
    slimit_max = n/10;
    slimit_min = n/40;
    slimit_step = (3.5*(slimit_max - slimit_min))/n;
   }
 else if (n >= 4000 && n < 8000)
   {
    slimit_max = n/25;
    slimit_min = n/125;
    slimit_step = (3.*(slimit_max - slimit_min))/n;
   }
 else if (n >= 8000 && n < 16000)
   {
    slimit_max = n/100;
    slimit_min = n/500;
    slimit_step = (2.5*(slimit_max - slimit_min))/n;
   }
 else if (n >= 16000 && n < 32000)
   {
    slimit_max = sqrt((REAL)n);
    slimit_min = n/1600;
    slimit_step = (2.0*(slimit_max - slimit_min))/n;
   }
 else if (n >= 32000 && n < 64000)
   {
    slimit_max = sqrt((REAL)n)/2.;
    slimit_min = 10;
    slimit_step = (1.75*(slimit_max - slimit_min))/n;
   }
 else
   {
    slimit_max = sqrt((REAL)n)/2.;
    slimit_min = 10;
    slimit_step = (1.5*(slimit_max - slimit_min))/n;
   }
/* args_fprintf(stdout, "slimits(%~,%~,%~)", slimit_min, slimit_max, slimit_step);*/
}

void raytrace(Screen* s)
{
 Ray r;
 Triangle* ts;
 int i,j,k;
 int idx;
 int prep_loaded;
 load_scene(&ts, "scene.4dd");
 vx_ndiv = 2;
 vx_whendiv = 5000;
 i = vx_whendiv;
 compute_slimits(i);
 
 btree = NULL;
 if (!(prep_loaded = load_preprocessed(ts, 1))) preprocess_scene(ts);
 free_mem();
 btree_info();
 g_ts = ts;
 if (!prep_loaded) save_preprocessed();
 loaded = 1;
 get_normal = get_normal_4d;
 intersection = intersection_4d;
 intersection_itree = intersection_itree_4d;
 r.P.x = observer.x;
 r.P.y = observer.y;
 r.P.z = observer.z;
 r.r = 0;
 r.d.z = (s->x+s->y)*lookz;
 args_fprintf(stdout, "\n");
 line_idx = 0;
 itree = alloc_itree();
 idx = 0;
 for (i=idx;i<s->x;i++)
   {
    r.d.x = i - s->x/2.;
    line_idx = i;
    for (j=0;j<s->y;j++)
      {
       r.d.y = j - s->y/2.;
       calculate_color(s, &r, ts, i, j);
      }
   }

/* wrt_bmp(&screen, screenbmp);*/
 free_scene(&ts);
}

int main(int lb, char** par)
{
 raytrace(&screen);
 free_screen(&screen);
 if (boxes) blist_free(&boxes);
 
 return 0;

}
