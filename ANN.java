import java.util.*;
import java.util.Scanner;
import java.io.*;
import java.io.File;
import java.io.FileWriter;

class Matrix {
  /* Base Matrix class borrowed from Math Department of Kentucky University
     Added function to matrix class in order to emulate popular R matrix capabilities such as 
     subsetting and dot.operations. */

  public int rows, columns;
  public double[][] element;  //the array containing the matrix

  public Matrix Copy() {
    Matrix result =  new Matrix(rows, columns);
    result.element = element;
    return result;
  }

  public Matrix() {
    // Create a zero 1x1 matrix
    rows = 1;
    columns = 1;
    element = new double[rows][columns];
  }

  public Matrix(int r, int c) {
    //  creates an empty r by c matrix
    rows = r;
    columns = c;
    element = new double[rows][columns];
  }

  public Matrix(double d) {
    //  creates a 1x1 matrix of double d
    rows = 1;
    columns = 1;
    element = new double[1][1];
    element[0][0] = d;
  }

  public Matrix(int r, int c, double fill) {
    //  creates an  r by c matrix with entries 'fill'
    rows = r;
    columns = c;
    element = new double[rows][columns];
    int i, j;
    for (i=0; i<rows; i++) {
      for (j=0; j<columns; j++) {
        element[i][j] = fill;
      }
    }
  }

  public Matrix(Matrix m) {
    //  creates a new replicate of m
    rows = m.rows;
    columns = m.columns;
    element = new double[rows][columns];
    int i, j;
    for (i=0; i<rows; i++) {
      for (j=0; j<columns; j++) {
        element[i][j] = m.element[i][j];
      }
    }
  }

  public Matrix(int r, int c, char code) {
  // contsructor: creates an  r by c special matrix
      rows = r;
      columns = c;
      element = new double[rows][columns];
      int i, j;

      if ((code == 'i') || (code == 'I')){
        // make an identity matrix
        for (i = 0; i < r; i++) {
            if (i < c) {
              element[i][i] = 1;
            }
        }
      }
      else if ((code == 'h') || (code == 'H')){
        // make a Hilbert matrix
        for (i = 0; i < r; i++) {
            for (j=0; j<c; j++) {
              element[i][j] = 1/((double)i+(double)j+1);
            }
        }
      }
      else if ((code == 'r') || (code == 'R')){
        // make a random matrix with entries uniform in [0, 1]
        for (i = 0; i < r; i++) {
            for (j=0; j<c; j++) {
              element[i][j] = Math.random()/1000;
            }
        }
      }
  }

  public Matrix transpose() { // transpose might not work
  // returns the transpose of this matrix object
      Matrix t = new Matrix(columns, rows);
      int i, j;
      for (i = 0; i<rows; i++) {
        for (j = 0; j<columns; j++) {
          t.element[j][i] = this.element[i][j];
        }
      }
      return t;
  }

  public static Matrix add(Matrix m1, Matrix m2){
    // Return the matrix m = m1 + m2
    Matrix m=new Matrix(m1.rows,m1.columns);
    if ((m1.rows == m2.rows)&&(m1.columns==m2.columns)) {
      int i,j;
      for (i=0; i<m.rows; i++) {
          for (j=0; j<m.columns; j++) {
            m.element[i][j] = m1.element[i][j] + m2.element[i][j]; 
          }
      }
    } else {
            System.out.println("Matrix dimension mismatch in add(m1,m2)");
            return new Matrix(1,1,1f);
        }
    return m;
  }

  public Matrix add(Matrix m2){
  //this + m2
      Matrix m = new Matrix(this.rows, this.columns);
      if ((this.rows == m2.rows) && (this.columns==m2.columns)) {
        int i,j;
        for(i = 0; i < m.rows; i++) {
            for(j = 0; j < m.columns; j++) {
              this.element[i][j] = this.element[i][j] + this.element[i][j]; 
            }
        }
      } else {
          System.out.println("Matrix dimension mismatch in m1.add(m2)");
          return new Matrix(1,1,1f);
      }
      return m;
  }

  public static Matrix subtract(Matrix m1, Matrix m2){
    // Return the matrix m = m1 - m2
    Matrix m=new Matrix(m1.rows,m1.columns);
    if ((m1.rows == m2.rows)&&(m1.columns==m2.columns)) {
      int i,j;
      for (i=0; i<m.rows; i++) {
        for (j=0; j<m.columns; j++) {
            m.element[i][j] = m1.element[i][j] - m2.element[i][j]; 
        }
      }
    } else {
        System.out.println("Matrix dimension mismatch subtract(m1, m2)");
        return new Matrix(1,1,1f);
    }
    return m;
  }

  public Matrix subtract(Matrix m2){
  // Return the matrix m = m1 - m2
      Matrix m=new Matrix(this.rows, this.columns);
      if ((this.rows == m2.rows)&&(this.columns==m2.columns)) {
        int i,j;
        for (i=0; i < m.rows; i++) {
            for (j=0; j < m.columns; j++) {
              m.element[i][j] = this.element[i][j] - m2.element[i][j]; 
            }
        }
      } else {
          System.out.println("Matrix dimension mismatch m1.subtract(m2)");
          return new Matrix(1,1,1f);
      }
      return m;
  }

  public static Matrix multiply(double d, Matrix m1){
    // Return the matrix m = d*m1
    Matrix m=new Matrix(m1.rows,m1.columns);
    int i,j;
    for (i=0; i<m.rows; i++) {
      for (j=0; j<m.columns; j++) {
        m.element[i][j] = d * m1.element[i][j];
      }
    }
    return m;
  }

  public static Matrix multiply(Matrix m1, Matrix m2){
      /* Matrix-Matrix or Matrix-vector product
         returns m=m1*m2
         m1 can be a 1x1 Matrix for scalar-Matrix product
      */
      Matrix m = new Matrix(0);
      if (m1.columns == m2.rows) {
        // matrix product
        double sum = 0;
        int k = 0;
        m = new Matrix(m1.rows,m2.columns);
        int i,j;
        for (i=0; i<m.rows; i++) {
            for (k=0; k<m2.columns; k++) {
              for (j=0; j<m1.columns; j++) {
                sum = sum + m1.element[i][j] * m2.element[j][k];
              }
              m.element[i][k] = sum;
              sum = 0;
            }
        }
      }
      else {
          System.out.println("dimension mismatch in multiply(m1, m2)");
      }
      return m;
  }

  public Matrix multiply(Matrix m2) {
      Matrix m = new Matrix(0);
      if (this.columns == m2.rows) {
        // matrix product
        double sum = 0;
        int k = 0;
        m = new Matrix(this.rows,m2.columns);
        int i,j;
        for (i=0; i<m.rows; i++) {
            for (k=0; k<m2.columns; k++) {
              for (j=0; j < this.columns; j++) {
                sum = sum + this.element[i][j] * m2.element[j][k];
              }
              m.element[i][k] = sum;
              sum = 0;
            }
        }
      }
      else {
          System.out.println("dimention mismatch in m1.multiply(m2)");
      }
      return m;
  }

  public Matrix subsetRows(int[] rows) {
      int nrows = rows.length;
      Matrix result = new Matrix(nrows, columns);
      for(int i = 0; i < nrows; ++i) {
          if(i >= this.rows) {
              System.out.println("invalid row number during subsetRows");
              return new Matrix(1,1,0);
          }
          result.element[i] = this.element[rows[i]];
      }
      return result;
  }

  public Matrix subsetRows(int start, int end) {
      int nrows = end-start+1;
      Matrix result = new Matrix(nrows, columns);
      for(int i = 0; i < nrows; ++i) {
          if(i >= this.rows) {
              System.out.println("invalid row number during subsetRows");
              return new Matrix(1,1,0);
          }
          result.element[i] = this.element[start+i];
      }
      return result;
  }

  public static Matrix dotMultiply(Matrix m1, Matrix m2) {
      if(m1.rows != m2.rows || m1.columns != m2.columns) {
        System.out.println("dimention mismatch in dotMultiply");
        return new Matrix(1,1,0);
      }
      Matrix result = new Matrix(m1.rows, m1.columns);
      for(int i = 0; i < m1.rows; ++i) {
        for(int j = 0; j < m1.columns; ++j) {
          result.element[i][j] = m1.element[i][j] * m2.element[i][j];
        }
      }
      return result;
  }

  public Matrix dotMultiply(Matrix m2) {
      if(this.rows != m2.rows || this.columns != m2.columns) {
        System.out.println("dimention mismatch in dotMultiply");
        return null;
      }
      Matrix result = new Matrix(this.rows, this.columns);
      for(int i = 0; i < this.rows; ++i) {
        for(int j = 0; j < this.columns; ++j) {
          result.element[i][j] = this.element[i][j] * m2.element[i][j];
        }
      }
      return result;
  }
  
  public Matrix appendCols(Matrix x){
    // append the column vectors in x to this matrix
    Matrix M = new Matrix(rows,columns+x.columns);
    int i,j;
    for(i=0;i<rows;i++){
      for(j=0;j<columns;j++){
        M.element[i][j]=this.element[i][j];
      }
      for(j=0;j<x.columns;j++){
        M.element[i][columns+j]=x.element[i][j];
      }
    }
    return M;
  }

  public Matrix prependCols(Matrix x) {
    // append the column vectors in x to this matrix
    Matrix M = new Matrix(rows,columns+x.columns);
    int i,j;
    for(i=0;i<rows;i++){
      for(j=0 ;j < x.columns; j++) {
        M.element[i][j] = x.element[i][j];
      }
      for(j=0; j < columns; j++) {
        M.element[i][x.columns+j] = this.element[i][j];
      }
    }
    return M;
  }

  public void Standardize() {
    double max = max();
    double min = min();
    double range = max - min;
    for(int i = 0; i < rows; ++i){
      for(int j = 0; j < columns; ++j){
        element[i][j] = element[i][j]/range;
      }
    }
  }
  
  public Matrix appendRows(Matrix x){
    // append the row vectors in x to this matrix
    Matrix M=new Matrix(rows+x.rows,columns);
    int i,j;
    for(i=0;i<columns;i++){
      for(j=0;j<rows;j++){
        M.element[j][i]=this.element[j][i];
      }
      for(j=0;j<x.rows;j++){
        M.element[rows+j][i]=x.element[j][i];
      }
    }
    return M;
  }

  public Matrix permute(int a1, int a2, char c) {
    /*  Returns a permuted matrix according  code c
    where c is in {'c', 'r'} for columns or rows and
    a1, a2 represent the columns/rows to swap
    */
    Matrix p = new Matrix(this);
    int i, j;
    if (c == 'r') {
      for (i=0; i<columns; i++) {
        p.element[a1][i] = this.element[a2][i];
        p.element[a2][i] = this.element[a1][i];
      }
    }
    else if (c == 'c') { 
      for (i=0; i<rows; i++) {
        p.element[i][a1] = this.element[i][a2];
        p.element[i][a2] = this.element[i][a1];
      }
    }
    return p;
  }

  public double norm() {
    /* returns the Frobenius norm (Matrix), or Euclidean norm (Vector)
       This is the default norm for a Matrix object. Use the Norm
       class for different norms.
    */
    double l = 0;
    int i, j;
    for (i = 0; i<rows; i++) {
      for (j = 0; j<columns; j++) {
        l = l + this.element[i][j] * this.element[i][j];
      }
    }
    l = Math.pow(l, 0.5);
    return l;
  }
  
  public double max() {
    /* returns the most positive element of the matrix or vector */
    double m = this.element[0][0];
    int i,j;
    for (i = 0; i<rows; i++) {
      for (j = 0; j<columns; j++) {
        if(this.element[i][j] > m){
          m = this.element[i][j];
        }
      }
    }
    return m;
  }
  
  public double min() {
    /* returns the most positive element of the matrix or vector */
    double m = this.element[0][0];
    int i,j;
    for (i = 0; i<rows; i++) {
      for (j = 0; j<columns; j++) {
        if(this.element[i][j] < m){
          m = this.element[i][j];
        }
      }
    }
    return m;
  }

  public double sum() {
    /* returns the sum of all the elements in the matrix or vector
     */
    double s = 0;
    int i, j;
    for (i = 0; i<rows; i++) {
      for (j = 0; j<columns; j++) {
        s = s + this.element[i][j];
      }
    }
    return s;
  }
  
  public double average() {
    // returns the average of all the elements in the matrix or vector
    double s = 0;
    int i, j;
    for (i = 0; i<rows; i++) {
      for (j = 0; j<columns; j++) {
        s = s + this.element[i][j];
      }
    }
    return s/(columns*rows);
  }
  
  public double sumSquares() {
    // returns the sum of the squares 
    // of all the elements in the matrix or vector
    double s = 0;
    int i, j;
    for (i = 0; i<rows; i++) {
      for (j = 0; j<columns; j++) {
        s = s + Math.pow(this.element[i][j],2);
      }
    }
    return s;
  }

  public String getElements() {
    String result = "";
    for(int i = 0; i < rows; ++i) {
      for(int j = 0; j < columns; ++j) {
        result += element[i][j] + " ";
      }
      result += "\n";
    }
    return result;
  }
}


public class ANN{

    public ArrayList<Matrix> layers = new ArrayList<Matrix>();
    public ArrayList<Matrix> aValues = new ArrayList();
    public ArrayList<Matrix> zValues = new ArrayList();
    public static Map<Matrix,Matrix> mappingXnY = new HashMap<Matrix,Matrix>();
    public static ArrayList<Matrix> matrixList = new ArrayList<Matrix>();

    public static void println(String s) {
      System.out.println(s);
    }

    public ANN(int[] layerSizes){
        for(int i = 0; i < layerSizes.length -1; ++i) {
            Matrix curTheta;
            curTheta = new Matrix(1 + layerSizes[i], layerSizes[i+1], 'r');
            layers.add(curTheta);
        }
    }

    public static Matrix parseFile(String filename) {
      try {
        Scanner parser = new Scanner(new File(filename));
        Matrix result = new Matrix(1, 128 * 120);
        int counter = 0;
        while(parser.hasNext()) {
          result.element[0][counter] = parser.nextDouble();
          ++counter;
        }
        if(counter != 128 * 120) {
          println("error during input");
          return null;
        } else {
          return result;
        }
      } catch (FileNotFoundException e) {
        println("file not found: " + filename);
        return null;
      }
    }

    public Matrix Predict(Matrix input) {
       aValues.clear();
       zValues.clear();

       int m = input.rows;
       Matrix one = new Matrix(m,1,1.0);

       Matrix curLayer = input.Copy();
       curLayer = curLayer.prependCols(one);
       aValues.add(curLayer); 

       curLayer = curLayer.multiply(layers.get(0));
       zValues.add(curLayer);
       curLayer = sigmoid(curLayer);
       curLayer = curLayer.prependCols(one);
       aValues.add(curLayer);
 

       for(int i = 1; i < layers.size(); ++i) {
           curLayer = Matrix.multiply(curLayer, layers.get(i));
           zValues.add(curLayer);
           curLayer = sigmoid(curLayer);
           curLayer = curLayer.prependCols(one);
           aValues.add(curLayer);
        }
       return curLayer;
    }

    public Matrix PredictFilter(Matrix input) {
      input = input.transpose();
      input = input.subsetRows(1, 2);
      input = input.transpose();
      for(int i = 0; i < input.rows; ++i) {
        if(input.element[i][0] > input.element[i][1]) {
          input.element[i][0] = 1;
          input.element[i][1] = 0;
        } else {
          input.element[i][0] = 0;
          input.element[i][1] = 1;
        }
      }
      return input;
    }

    public Matrix sigmoid(Matrix input) {
      for(int i = 0; i < input.rows; ++i) {
          for(int j = 0; j < input.columns; ++j) {
              input.element[i][j] = 1 / (1 + Math.exp(input.element[i][j] * -1));
          }
      }
        return input; // check here
    }

    public Matrix sigmoidGradient(Matrix input) {
      Matrix sigmoid = sigmoid(input);
      Matrix ones = new Matrix(input.rows, input.columns, 1.0);
      Matrix gradient = Matrix.dotMultiply(sigmoid, Matrix.subtract(ones, sigmoid));
      return gradient;
    }

    public void dim(Matrix m) {
      println(m.rows + " " + m.columns + " ");
    }

    public void BP(Matrix pred, Matrix Y, float alpha) {
      Matrix predTemp = pred;
      predTemp = predTemp.transpose();
      predTemp = predTemp.subsetRows(1, predTemp.rows-1);

      Matrix outputDelta = Matrix.subtract(predTemp, Y.transpose());
      int m = Y.rows;

      ArrayList<Matrix> gradients = new ArrayList<Matrix>();
      ArrayList<Matrix> deltas = new ArrayList<Matrix>();

      deltas.add(outputDelta);

      for(int i = layers.size()-1; i >= 1; --i) { // each layer
        Matrix curDeltatemp = layers.get(i).multiply(deltas.get(deltas.size()-1));
        curDeltatemp = curDeltatemp.subsetRows(1, curDeltatemp.rows-1);

        Matrix curDelta = curDeltatemp.dotMultiply(sigmoidGradient(zValues.get(i-1)).transpose());

        deltas.add(curDelta);
      }

      Collections.reverse(deltas);
      for(int i = 0; i < deltas.size(); ++i) {
        Matrix curGradient = Matrix.multiply((float)1/m * alpha, aValues.get(i).transpose().multiply(deltas.get(i).transpose()));
        if(i == 1) {
        }
        gradients.add(curGradient);
      }

      for(int i = 0; i < gradients.size(); ++i) {
        layers.set(i, layers.get(i).subtract(gradients.get(i)));
      }
    }

    public static void main(String[] args) {
      int[] layerSizes = {128 * 120, 200, 2};
      ANN neuralNet = new ANN(layerSizes);

      Matrix bigX = new Matrix(); 
    Matrix femaleMatrix = new Matrix(1,2,0);
    femaleMatrix.element[0][1]=1; 
    
    Matrix maleMatrix = new Matrix(1,2,0);
    maleMatrix.element[0][0]=1; 

    Matrix bigY = femaleMatrix;
    
    String fulldir = (System.getProperty("user.dir") + "/Female/");
    File folder = new File(fulldir);
    String[] filenames = folder.list(new FilenameFilter() {
      @Override
          public boolean accept(File dir, String name) {
              return !name.equals(".DS_Store");
          }
      });
    
    bigX = ANN.parseFile(fulldir + filenames[0]);
    mappingXnY.put(bigX, femaleMatrix);

    for(int j = 1; j < filenames.length; j++)
    {
      Matrix temp = neuralNet.parseFile(fulldir + filenames[j]);
      bigX = bigX.appendRows(temp);
      mappingXnY.put(temp, femaleMatrix);
      bigY = bigY.appendRows(femaleMatrix); 
    }

    fulldir = System.getProperty("user.dir")+"/Male/";
    folder = new File(fulldir);

    filenames = folder.list(new FilenameFilter() {
      @Override
          public boolean accept(File dir, String name) {
              return !name.equals(".DS_Store");
          }
      });

    for(int j = 0; j < filenames.length; j++)
    {
      Matrix temp = ANN.parseFile(fulldir + filenames[j]);
      bigX=bigX.appendRows(temp);
      mappingXnY.put(temp,maleMatrix);
      bigY=bigY.appendRows(maleMatrix);
    }

    try
    {
      int i = 0;
      while(i < args.length)
      {
        if(args[i].equalsIgnoreCase("-train"))
        {
          println("in train");

            println(bigY.getElements());
            
            println("begin training");
            neuralNet.dim(bigX);
            Matrix pred = neuralNet.Predict(bigX);

            int trainings = 10;
            while(trainings > 0) {
              println("training # " + trainings);
              neuralNet.BP(pred,bigY,3f);
              pred = neuralNet.Predict(bigX);
              trainings--;
            }
            println(pred.getElements());
            pred = neuralNet.PredictFilter(pred);
            println(pred.getElements());

            Matrix diff = pred.subtract(bigY);
            println(diff.getElements());
        }

        else if(args[i].equalsIgnoreCase("-cv")) {
          neuralNet.crossValidation();
        }
        else if(args[i].equalsIgnoreCase("-test"))
        {
          println("begin training");
            neuralNet.dim(bigX);
            Matrix pred = neuralNet.Predict(bigX);

            int trainings = 25;
            while(trainings > 0) {
              println("training # " + trainings);
              neuralNet.BP(pred,bigY,3f);
              pred = neuralNet.Predict(bigX);
              trainings--;
            }
            println(pred.getElements());
            pred = neuralNet.PredictFilter(pred);
            println(pred.getElements());

            Matrix diff = pred.subtract(bigY);
            println(diff.getElements());

            /*testing */
            bigX = ANN.parseFile(fulldir + filenames[0]);
          for(int j = 1; j < filenames.length; j++)
          {
            Matrix temp = ANN.parseFile(fulldir + filenames[j]);
            bigX = bigX.appendRows(temp);
          }
          pred = neuralNet.Predict(bigX);
            try{
            FileWriter fr = new FileWriter("998241764_999376654.prediction");
            BufferedWriter br = new BufferedWriter(fr);
            PrintWriter out = new PrintWriter(br);
            for(int r = 0; r < pred.rows; ++r) {
              double c = Math.abs(pred.element[r][1]-pred.element[r][2])/1.0;
              out.write("{"+ pred.element[r][1] +","+ pred.element[r][2] +"}"+" "+c);
              out.write(" ");
              out.write("\n");
            }
            out.close();
          }catch(IOException e)
          {
            System.out.println(e);
          }

          /*heat map */
          try{
            FileWriter fr = new FileWriter("weightmap.txt");
            BufferedWriter br = new BufferedWriter(fr);
            PrintWriter out = new PrintWriter(br);
            for(int r = 0; r < neuralNet.layers.get(i).rows; ++r) {
                out.write(Double.valueOf(neuralNet.layers.get(i).element[r][0]).toString() + " ");
            }
            out.close();
          }catch(IOException e){
            System.out.println(e);
          }
            
            fulldir = System.getProperty("user.dir")+"/Test/";
          folder = new File(fulldir);
            filenames = folder.list(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return !name.equals(".DS_Store");
                }
            });
            for(int j = 0; j < filenames.length; j++)
            {
              Matrix temp = ANN.parseFile(fulldir+filenames[j]);
              ANN.println("in loop");
              ANN.println(filenames[j]);
            }
        }
        else 
          throw new IllegalArgumentException();
        i += 2;
      }
    }
    catch(IndexOutOfBoundsException ioob)
    {
      System.err.println("Invalid Arguments");
      System.exit(2);
    }
    catch(IllegalArgumentException ia)
    {
      System.err.println("Invalid Arguments: " + ia.getMessage());
      System.exit(4);
    }
    catch(Exception e)
    {
      System.err.println("Unknown Error");
      System.exit(5);
    }
  }
    
  public static void crossValidation() {
    ArrayList<Matrix> xList = new ArrayList<Matrix>(mappingXnY.keySet());
    List<Matrix> testList = new ArrayList<Matrix>();
    List<Matrix> trainList;
    Matrix bigX;
    Matrix bigY;  

    //shuffle the original data
    Collections.shuffle(xList);

    System.out.println("data size"+mappingXnY.size());
    int dataSize = mappingXnY.size();
    int setSize = dataSize/5;
    int extras = dataSize%5;

    int numExperiment = 30;
    double[] trainingAccuracyExp = new double[5];
    double[] testingAccuracyExp = new double[5];
    double accuracy = 0;

    while(numExperiment > 0) {
      numExperiment--; 
      for(int i = 0; i < 5; i++){
        bigX = new Matrix();
        bigY = new Matrix();
        trainList = new ArrayList<Matrix>(xList);
        if(i < extras){
          testList.addAll(trainList.subList(i*(setSize+1),(i+1)*(setSize+1)));
          trainList.subList(i*(setSize+1),(i+1)*(setSize+1)).clear();
          System.out.println("num:"+i);
          System.out.println("size of test:"+testList.size());
          System.out.println("size of train:"+trainList.size());    
        }
        else{
          testList.addAll(trainList.subList(i*setSize,(i+1)*setSize));
          trainList.subList(i*setSize,(i+1)*setSize).clear();
          System.out.println("num:"+i);
          System.out.println("size of test:"+testList.size());
          System.out.println("size of train:"+trainList.size());
        }

        bigX = trainList.get(0);
        bigY = mappingXnY.get(bigX);

        for(int j = 1; j < trainList.size();j++){
          bigX = bigX.appendRows(trainList.get(j));
          bigY = bigY.appendRows(mappingXnY.get(trainList.get(j)));
        }

        int[] layerSizes = {128 * 120, 200, 2};
        ANN neuralNet = new ANN(layerSizes);
        neuralNet.dim(bigX);
        neuralNet.dim(bigY);

        println("begin training");
        Matrix pred = neuralNet.Predict(bigX);

        int trainings = 25;
        while(trainings > 0) {
          println("training # " +trainings);
          neuralNet.BP(pred,bigY,3f);
          pred = neuralNet.Predict(bigX);
          trainings--;
        }

        pred = neuralNet.PredictFilter(pred);
        Matrix diff = bigY.subtract(pred);
        int correct = 0;
        for(int j = 0; j < diff.rows; ++j) {
          if(diff.element[j][0] == 0) {
            correct++;
          }
        }
        System.out.println("corrects: " + correct);
        System.out.println("accuracy = " + (float)((float)correct/(float)diff.rows));
        accuracy = (float)((float)correct/(float)diff.rows);
        trainingAccuracyExp[i] = accuracy;
        
        accuracy = 0;
        bigX = new Matrix();
        bigY = new Matrix();
        bigX = testList.get(0);
        bigY = mappingXnY.get(bigX);
        for(int j = 1; j < testList.size();j++){
          bigX = bigX.appendRows(testList.get(j));
          bigY = bigY.appendRows(mappingXnY.get(testList.get(j)));
        }


        pred = neuralNet.Predict(bigX);
        pred = neuralNet.PredictFilter(pred);
        diff = bigY.subtract(pred);
        correct = 0;
        for(int j = 0; j < diff.rows; ++j) {
          if(diff.element[j][0] == 0) {
            correct++;
          }
        }
        System.out.println("corrects: " + correct);
        System.out.println("accuracy = " + (float)((float)correct/(float)diff.rows));
        accuracy = (float)((float)correct/(float)diff.rows);

        testingAccuracyExp[i] = accuracy;
        testList.clear();
        trainList.clear();
      }
      double mean = 0;
      double sd = 0;
      for(int i = 0; i < 5; ++i) {
        mean += trainingAccuracyExp[i];
      }

      mean = (double)(mean / 5.0);

      for(int i = 0; i < 5; ++i) {
        sd += (trainingAccuracyExp[i]-mean) * (trainingAccuracyExp[i]-mean);
      }

      sd = Math.sqrt(sd);
      System.out.println("experiment #" + (10 - numExperiment));
      System.out.println("training mean " + mean);
      System.out.println("training sd : " + sd);

      mean = 0;
      sd = 0;
        
      for(int i = 0; i < 5; ++i) {
        mean += trainingAccuracyExp[i];
      }

      mean = mean / 5.0;

      for(int i = 0; i < 5; ++i) {
        sd += (testingAccuracyExp[i]-mean) * (testingAccuracyExp[i]-mean);
      }
      sd = Math.sqrt(sd);
      System.out.println("testing mean " + mean);
      System.out.println("testing sd : " + sd);

    }
    }
}
