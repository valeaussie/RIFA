#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "shp3.h"

int read_file_header(FILE *in, ShpHeader *head)
{
	int i;
	char *pch;

	pch = (char*)(&(head->file_code));
	for (i = 0; i < 4; i++) {
		fscanf(in, "%c", pch + (3 - i));
	}

	for (i = 0; i < 20; i++) {
		fscanf(in, "%c", head->blank + i);
	}

	pch = (char*)(&(head->length));
	for (i = 0; i < 4; i++) {
		fscanf(in, "%c", pch + (3 - i));
	}

	pch = (char*)(&(head->version));
	for (i = 0; i < 4; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->shape_type));
	for (i = 0; i < 4; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->x_min));
	for (i = 0; i < 8; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->y_min));
	for (i = 0; i < 8; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->x_max));
	for (i = 0; i < 8; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->y_max));
	for (i = 0; i < 8; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->z_min));
	for (i = 0; i < 8; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->z_max));
	for (i = 0; i < 8; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->m_min));
	for (i = 0; i < 8; i++) {
		fscanf(in, "%c", pch + i);
	}

	pch = (char*)(&(head->m_max));
	for (i = 0; i < 8; i++) {
		fscanf(in, "%c", pch + i);
	}

	if (head->shape_type != 1 && head->shape_type != 5) {
		printf("Shape file must contain point type shapes (1) or polygon type shapes (5).\n");
		exit(1);
	}

	return 0;
	//TO DO: return 1 if any of these reads fail
}


int output_file_header(FILE *out, ShpHeader *head)
{
	int i;
	char *pch;

	pch = (char*)(&(head->file_code));
	for (i = 0; i < 4; i++) {
		fprintf(out, "%c", pch[3 - i]);
	}

	for (i = 0; i < 20; i++) {
		fprintf(out, "%c", head->blank[i]);
	}

	pch = (char*)(&(head->length));
	for (i = 0; i < 4; i++) {
		fprintf(out, "%c", pch[3 - i]);
	}

	pch = (char*)(&(head->version));
	for (i = 0; i < 4; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->shape_type));
	for (i = 0; i < 4; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->x_min));
	for (i = 0; i < 8; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->y_min));
	for (i = 0; i < 8; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->x_max));
	for (i = 0; i < 8; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->y_max));
	for (i = 0; i < 8; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->z_min));
	for (i = 0; i < 8; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->z_max));
	for (i = 0; i < 8; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->m_min));
	for (i = 0; i < 8; i++) {
		fprintf(out, "%c", pch[i]);
	}

	pch = (char*)(&(head->m_max));
	for (i = 0; i < 8; i++) {
		fprintf(out, "%c", pch[i]);
	}

	return 0;
	//TO DO: return 1 if any of these reads fail
}


int read_record_header(FILE *in, ShpRecHeader *rechead)
{
	int i;
	char *pch;

	pch = (char*)(&(rechead->record_number));
	for (i = 0; i < 4; i++) {
		if (fscanf(in, "%c", pch + (3 - i)) != 1) return 1;
	}

	pch = (char*)(&(rechead->content_length));
	for (i = 0; i < 4; i++) {
		if (fscanf(in, "%c", pch + (3 - i)) != 1) return 1;
	}

	return 0;
}


int output_record_header(FILE *out, ShpRecHeader *rechead)
{
	int i;
	char *pch;

	pch = (char*)(&(rechead->record_number));
	for (i = 0; i < 4; i++) {
		if (fprintf(out, "%c", pch[3 - i]) != 1) return 1;
	}

	pch = (char*)(&(rechead->content_length));
	for (i = 0; i < 4; i++) {
		if (fprintf(out, "%c", pch[3 - i]) != 1) return 1;
	}

	return 0;
}


int read_point_record(FILE *in, ShpPointRecord *rec)
{
	int i;
	char *pch;

	pch = (char*)(&(rec->shape_type));
	for (i = 0; i < 4; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	if (rec->shape_type == 0) return 0; //Null shape

	pch = (char*)(&(rec->x));
	for (i = 0; i < 8; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	pch = (char*)(&(rec->y));
	for (i = 0; i < 8; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	return 0;
}

void read_points(char *filename, ShpHeader *head, ShpRecHeader **recheads, ShpPointRecord **records, int *numrec)
{
	FILE *in;
	ShpRecHeader testrechead;
	ShpPointRecord testrec;
	int i;

	in = fopen(filename, "r");
	if (!in) { printf("Could not open file %s.\n", filename); exit(1); }

	read_file_header(in, head);

	//Count records
	*numrec = 0;
	while (read_record_header(in, &testrechead) == 0) {
		read_point_record(in, &testrec); //Assumes it won't fail!
		(*numrec)++;
	}
	printf("Number of point records is %d.\n", *numrec);

	//Read records
	rewind(in);
	read_file_header(in, head);
	*recheads = (ShpRecHeader*)malloc((*numrec) * sizeof(ShpRecHeader));
	if (!(*recheads)) { printf("malloc error in read_points.\n"); exit(1); }
	*records = (ShpPointRecord*)malloc((*numrec) * sizeof(ShpPointRecord));
	if (!(*records)) { printf("malloc error in read_points.\n"); exit(1); }
	for (i = 0; i < *numrec; i++) {
		read_record_header(in, *recheads + i);
		read_point_record(in, *records + i);
	}

	fclose(in);
}


void cleanup_points(int numrec, ShpRecHeader *recheads, ShpPointRecord *records)
{
	free(recheads);
	free(records);
}


int read_poly_record(FILE *in, ShpPolygonRecord *rec)
{
	int i, j;
	char *pch;

	pch = (char*)(&(rec->shape_type));
	for (i = 0; i < 4; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	if (rec->shape_type == 0) return 0; //Null shape

	pch = (char*)(&(rec->x_min));
	for (i = 0; i < 8; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	pch = (char*)(&(rec->y_min));
	for (i = 0; i < 8; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	pch = (char*)(&(rec->x_max));
	for (i = 0; i < 8; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	pch = (char*)(&(rec->y_max));
	for (i = 0; i < 8; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	pch = (char*)(&(rec->num_parts));
	for (i = 0; i < 4; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	pch = (char*)(&(rec->num_points));
	for (i = 0; i < 4; i++) {
		if (fscanf(in, "%c", pch + i) != 1) return 1;
	}

	rec->parts = (int*)malloc(rec->num_parts * sizeof(int));
	if (!rec->parts) { printf("malloc error in read_poly_record.\n"); exit(1); }
	for (j = 0; j < rec->num_parts; j++) {
		pch = (char*)(rec->parts + j);
		for (i = 0; i < 4; i++) {
			if (fscanf(in, "%c", pch + i) != 1) return 1;
		}
	}

	rec->points = (Point*)malloc(rec->num_points * sizeof(Point));
	if (!rec->points) { printf("malloc error in read_poly_record.\n"); exit(1); }
	for (j = 0; j < rec->num_points; j++) {
		pch = (char*)(rec->points + j);
		for (i = 0; i < 16; i++) {
			if (fscanf(in, "%c", pch + i) != 1) return 1;
		}
	}

	return 0;
}


int output_poly_record(FILE *out, ShpPolygonRecord *rec)
{
	int i, j;
	char *pch;

	pch = (char*)(&(rec->shape_type));
	for (i = 0; i < 4; i++) {
		if (fprintf(out, "%c", pch[i]) != 1) return 1;
	}

	if (rec->shape_type == 0) return 0; //Null shape

	pch = (char*)(&(rec->x_min));
	for (i = 0; i < 8; i++) {
		if (fprintf(out, "%c", pch[i]) != 1) return 1;
	}

	pch = (char*)(&(rec->y_min));
	for (i = 0; i < 8; i++) {
		if (fprintf(out, "%c", pch[i]) != 1) return 1;
	}

	pch = (char*)(&(rec->x_max));
	for (i = 0; i < 8; i++) {
		if (fprintf(out, "%c", pch[i]) != 1) return 1;
	}

	pch = (char*)(&(rec->y_max));
	for (i = 0; i < 8; i++) {
		if (fprintf(out, "%c", pch[i]) != 1) return 1;
	}

	pch = (char*)(&(rec->num_parts));
	for (i = 0; i < 4; i++) {
		if (fprintf(out, "%c", pch[i]) != 1) return 1;
	}

	pch = (char*)(&(rec->num_points));
	for (i = 0; i < 4; i++) {
		if (fprintf(out, "%c", pch[i]) != 1) return 1;
	}

	rec->parts = (int*)malloc(rec->num_parts * sizeof(int));
	if (!rec->parts) { printf("malloc error in read_poly_record.\n"); exit(1); }
	for (j = 0; j < rec->num_parts; j++) {
		pch = (char*)(rec->parts + j);
		for (i = 0; i < 4; i++) {
			if (fprintf(out, "%c", pch[i]) != 1) return 1;
		}
	}

	rec->points = (Point*)malloc(rec->num_points * sizeof(Point));
	if (!rec->points) { printf("malloc error in read_poly_record.\n"); exit(1); }
	for (j = 0; j < rec->num_points; j++) {
		pch = (char*)(rec->points + j);
		for (i = 0; i < 16; i++) {
			if (fprintf(out, "%c", pch[i]) != 1) return 1;
		}
	}

	return 0;
}


int bndcmp(const void *arg1, const void *arg2)
{
	band *b1, *b2;

	b1 = (band*)arg1;
	b2 = (band*)arg2;

	if (b1->min < b2->min) return -1;
	if (b1->min > b2->min) return 1;
	return 0;
}


void createbands(ShpPolygonRecord *poly)
{
	int i, j;
	band *ybands, *yb;
	int upper;
	double y;
	int a, b, c;

	ybands = (band*)malloc((poly->num_points) * sizeof(band));
	if (!ybands) { printf("malloc error in createbands.\n"); exit(1); }

	//Copy and sort
	for (i = 0; i < poly->num_points; i++)
		ybands[i].min = poly->points[i].y_coordinate;
	qsort(ybands, poly->num_points, sizeof(band), &bndcmp);

	//Remove duplicates and count
	j = 1;
	for (i = 1; i < poly->num_points; i++) {
		if (ybands[i].min > ybands[i - 1].min) {
			ybands[j++].min = ybands[i].min;
		}
	}
	poly->num_bands = j; //Note the maximum has its own band

	//Zero counts
	for (i = 0; i < poly->num_bands; i++)
		ybands[i].numseg = 0;

	//Scan segments counting the number crossing each band
	for (i = 0; i < poly->num_parts; i++) {
		//Find first point
		j = poly->parts[i];
		y = poly->points[j].y_coordinate;

		//binary search
		a = 0;
		b = poly->num_bands - 1;
		if (ybands[b].min <= y) a = b;
		else while (b > a + 1) {
			c = (a + b) / 2;
			if (ybands[c].min > y) b = c;
			else a = c;
		}
		yb = ybands + a;

		//Scan segments
		if (i == poly->num_parts - 1) upper = poly->num_points - 1;
		else upper = poly->parts[i + 1] - 1;
		for (j++; j <= upper; j++) {
			if (poly->points[j].y_coordinate > poly->points[j - 1].y_coordinate) {
				do {
					yb->numseg++;
					yb++;
				} while (yb->min < poly->points[j].y_coordinate);
			}
			else if (poly->points[j].y_coordinate < poly->points[j - 1].y_coordinate) {
				do {
					yb--;
					yb->numseg++;
				} while (yb->min > poly->points[j].y_coordinate);
			}
			else { //Horizontal line
				yb->numseg++;
			}
		}
	}

	//Allocate space for segment list
	for (i = 0; i < poly->num_bands; i++) {
		yb = ybands + i;
		yb->seg = (int *)malloc(yb->numseg * sizeof(int));
		if (!yb->seg) { printf("malloc error in createbands.\n"); exit(1); }
		yb->numseg = 0;
	}

	//Redo and record segments
	for (i = 0; i < poly->num_parts; i++) {
		//Find first point
		j = poly->parts[i];
		y = poly->points[j].y_coordinate;

		//binary search
		a = 0;
		b = poly->num_bands - 1;
		if (ybands[b].min <= y) a = b;
		else while (b > a + 1) {
			c = (a + b) / 2;
			if (ybands[c].min > y) b = c;
			else a = c;
		}
		yb = ybands + a;

		//Scan segments
		if (i == poly->num_parts - 1) upper = poly->num_points - 1;
		else upper = poly->parts[i + 1] - 1;
		for (j++; j <= upper; j++) {
			if (poly->points[j].y_coordinate > poly->points[j - 1].y_coordinate) {
				do {
					yb->seg[yb->numseg] = j - 1;
					yb->numseg++;
					yb++; //Note there is an extra band at the end holding y_max
				} while (yb->min < poly->points[j].y_coordinate);
			}
			else if (poly->points[j].y_coordinate < poly->points[j - 1].y_coordinate) {
				do {
					yb--;
					yb->seg[yb->numseg] = j - 1;
					yb->numseg++;
				} while (yb->min > poly->points[j].y_coordinate);
			}
			else {
				yb->seg[yb->numseg] = j - 1;
				yb->numseg++;
			}
		}
	}

	poly->ybands = ybands;
}


void read_polygons(char *filename, ShpHeader *head, ShpRecHeader **recheads, ShpPolygonRecord **records, int *numrec)
{
	FILE *in;
	ShpRecHeader testrechead;
	ShpPolygonRecord testrec;
	int i;

	in = fopen(filename, "r");
	if (!in) { printf("Could not open file %s.\n", filename); exit(1); }

	read_file_header(in, head);

	//Count records
	*numrec = 0;
	while (read_record_header(in, &testrechead) == 0) {
		if (read_poly_record(in, &testrec) == 1) {
			printf("Failure to read polygon record in read_polygons.\n");
			exit(1);
		}
		(*numrec)++;

		if (testrec.shape_type == 5) {
			free(testrec.parts);
			free(testrec.points);
		}
	}

	printf("Number of polygon records in %s is %d.\n", filename, *numrec);

	//Read records
	rewind(in);
	read_file_header(in, head);
	*recheads = (ShpRecHeader*)malloc((*numrec) * sizeof(ShpRecHeader));
	if (!(*recheads)) { printf("malloc error in read_polygons.\n"); exit(1); }
	*records = (ShpPolygonRecord*)malloc((*numrec) * sizeof(ShpPolygonRecord));
	if (!(*records)) { printf("malloc error in read_polygons.\n"); exit(1); }
	for (i = 0; i < *numrec; i++) {
		read_record_header(in, *recheads + i);
		read_poly_record(in, *records + i);
		if ((*records)[i].shape_type == 5) {
			createbands(*records + i);
		}
	}

	fclose(in);
}


void cleanup_polygons(int numrec, ShpRecHeader *recheads, ShpPolygonRecord *records)
{
	int i, j;

	for (i = 0; i < numrec; i++) {
		if (records[i].shape_type == 5) {
			free(records[i].parts);
			free(records[i].points);
			for (j = 0; j < records[i].num_bands; j++)
				free(records[i].ybands[j].seg);
			free(records[i].ybands);
		}
	}
	free(recheads);
	free(records);
}


int in_interior(Point *p, ShpPolygonRecord *poly)
//determines whether point is in interior of polygon - returns 1 if so, 0 otherwise
//Note that if p is on a band min one should probably check the segments from the band below also.
//For example, a single point on a band min will not be identified as interior. Does it matter?
{
	int i;
	Point *p1, *p2;
	double crossx;
	int a, b, c;
	band *yband;
	int leftcnt;

	if (poly->shape_type == 0) return 0; //Null shape

	if (p->x_coordinate < poly->x_min || p->x_coordinate > poly->x_max ||
		p->y_coordinate < poly->y_min || p->y_coordinate > poly->y_max)
		return 0; //Outside bounding box of polygon

	  //Binary search to find appropriate band - assumes polygon bounds are tight! Check this when creating bands.
	a = 0;
	b = poly->num_bands - 1;
	if (poly->ybands[b].min <= p->y_coordinate) a = b;
	else while (b > a + 1) {
		c = (a + b) / 2;
		if (poly->ybands[c].min > p->y_coordinate) b = c;
		else a = c;
	}
	yband = poly->ybands + a;

	//Count segments to left of p
	leftcnt = 0;
	for (i = 0; i < yband->numseg; i++) {
		p1 = poly->points + yband->seg[i];
		p2 = p1 + 1;

		//Treat horizontal lines separately
		if (p1->y_coordinate == p2->y_coordinate) { //Horizontal line
			if ((p->y_coordinate == p1->y_coordinate) &&
				((p1->x_coordinate <= p->x_coordinate && p->x_coordinate <= p2->x_coordinate) ||
				(p1->x_coordinate >= p->x_coordinate && p->x_coordinate >= p2->x_coordinate)))
				return 1; //Point is on segment
		} //Ignore horizontal lines if point is not on segment
		else {
			crossx = p1->x_coordinate + 1.0*(p2->x_coordinate - p1->x_coordinate)
				*(p->y_coordinate - p1->y_coordinate) / (p2->y_coordinate - p1->y_coordinate);
			if (crossx <= p->x_coordinate) leftcnt++;
		}
	}

	//Determine whether inside or outside
	if ((leftcnt / 2) * 2 != leftcnt) return 1;
	return 0;
}


int get_id(Point *p, ShpPolygonRecord *records, int numrec)
{
	int i;

	for (i = 0; i < numrec; i++) {
		if (in_interior(p, records + i)) return i;
	}

	return -1;
}

/*
 int main()
 {
	 char filename[] = "Positives.shp";
	 FILE *in,*out;
	 ShpHeader head;
	 ShpRecHeader *recheads;
	 ShpPointRecord *records;
	 int numrec,i;

	 in = fopen(filename,"r");
	 if(!in){printf("Could not open file %s.\n",filename);exit(1);}

	 read_file_header(in,&head);

	 read_points(filename,&head,&recheads,&records,&numrec);

	 fclose(in);

	 out = fopen("nest_coords_2020.csv","w");
	 if(!out){printf("Could not open file nest_coords_2020.csv.\n");exit(1);}

	 for(i=0;i<numrec;i++){
		 printf("%d,%lf,%lf\n",recheads[i].record_number,records[i].x,records[i].y);
		 fprintf(out,"%d,%lf,%lf\n",recheads[i].record_number,records[i].x,records[i].y);
	 }

	 fclose(out);

	 cleanup_points(numrec,recheads,records);

	 return 0;
 }

 int main()
 {
	 char filename[] = "Jobs_L.shp";
	 FILE *in,*out;
	 ShpHeader head;
	 ShpRecHeader *recheads;
	 ShpPolygonRecord *records;
	 int numrec,i;


	 read_polygons(filename,&head,&recheads,&records,&numrec);

	 cleanup_polygons(numrec,recheads,records);


	 return 0;
 }
 */
