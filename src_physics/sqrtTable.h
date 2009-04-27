//calculates divided difference to use for newton's method (yVals gets destroyed)
void divided_diff(BigReal *xVals, BigReal *yVals, BigReal *dd, int len){
  dd[0] = yVals[0];
  for (int i = 1; i < len; i++){
    for (int j = 0; j < len-i; j++){
      yVals[j] = (yVals[j+1]-yVals[j])/(xVals[j+i]-xVals[j]);
    }
    dd[i] = yVals[0];
  }
}

//get constants for Horner's method
void getHornerForm(BigReal *dd, BigReal *xVals, sqrtPars *ps){
  ps->d = dd[3];
  ps->c = -dd[3] * (xVals[0] + xVals[1] + xVals[2]) + dd[2];
  ps->b = dd[3] * (xVals[0] * (xVals[1] + xVals[2]) + xVals[1] * xVals[2]);
  ps->b += -dd[2] * (xVals[0] + xVals[1]) + dd[1];
  ps->a = dd[0] - dd[1] * xVals[0] + dd[2] * xVals[0] * xVals[1];
  ps->a -= dd[3] * xVals[0] * xVals[1] * xVals[2];
}

//interpolates at len segments using 3rd degree newtons method with evenly spaced intervals
sqrtTable* fillTable(int len, BigReal delt){
  BigReal prevX, prevY;
  BigReal *xVals;
  BigReal *yVals;
  BigReal *dd;
  sqrtTable *table = new (len) sqrtTable;
  table->length = len;
  table->delta = delt;

  xVals = new BigReal[4];
  yVals = new BigReal[4];
  dd = new BigReal[4];

  prevX = 0;
  prevY = 0;
  for (int i = 0; i < len; i++){
    xVals[0] = prevX;
    xVals[1] = prevX + delt / 3;
    xVals[2] = xVals[1] + delt / 3;
    xVals[3] = prevX + delt;
    prevX = prevX + delt;

    yVals[0] = prevY;
    for (int j = 1; j < 4; j++)
      yVals[j] = sqrt(xVals[j]);
    prevY = yVals[3];
    divided_diff(xVals, yVals, dd, 4);
    getHornerForm(dd, xVals, &(table->pars[i]));
  }    
  return table;
}


