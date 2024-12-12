#ifndef LORENTZ_HEADER
#define LORENTZ_HEADER

template <typename Field> struct LV {
  std::vector<Field> components = std::vector<Field>(4);

  LV(){
    components = std::vector<Field>(4);
  }

  LV(std::vector<Field> lv):components(lv){};

  LV(Field* lv) {
    for(int i = 0; i < 4; i++) components[i] = lv[i];
  }

  LV<Field> conj() const {
    LV<Field> output;
    for(int i = 0; i < 4; i++) output.components[i] = std::conj(components[i]);
    return output;
  }

  void print() const {
    std::cout << "{" << components[0] << ", "
                     << components[1] << ", "
                     << components[2] << ", "
                     << components[3] << "}" << std::endl;
  }

  // Define the multiplication operator for Field1 * LV<Field2>
  template <typename Field1>
  LV<typename std::common_type<Field1, Field>::type> operator*(const Field1& scalar) const {
    LV<typename std::common_type<Field1, Field>::type> result;
    for (int i = 0; i < 4; ++i) {
      result.components[i] = scalar*components[i];
    }
    return result;
  }

  Field* to_arr() const {
    return &components[0];
  }
};

template <typename Field> struct LM {
  std::vector<std::vector<Field>> components = std::vector<std::vector<Field>>(4, std::vector<Field>(4));
  LM(){
    components = std::vector<std::vector<Field>>(4, std::vector<Field>(4));
  }
  LM(std::vector<std::vector<Field>> lm):components(lm){};
  void print() const {
    std::cout << "{";
    for(int i = 0; i < 4; i++) {
      std::cout << "{";
      for(int j = 0; j < 4; j++) {
        std::cout << components[i][j] << ",";
      }
      std::cout << "}," << std::endl;
    }
    std::cout << std::endl;
  }

  // Define the multiplication operator for Field1 * LM<Field2>
  template <typename Field1>
  LM<typename std::common_type<Field1, Field>::type> operator*(const Field1& scalar) const {
    LM<typename std::common_type<Field1, Field>::type> result;
    for (int i = 0; i < 4; ++i) for(int j = 0; j < 4; j++) {
      result.components[i][j] *= scalar;
    }
    return result;
  }

  // Define the multiplication operator for LM<Field1> * LM<Field2>
  template <typename Field1>
  LM<typename std::common_type<Field1, Field>::type> operator*(const LM<Field1>& m) const {
    LM<typename std::common_type<Field1, Field>::type> result;
    for (int i = 0; i < 4; ++i) for(int j = 0; j < 4; j++) for(int k = 0; k < 4; k++) {
      result.components[i][j] += m.components[i][k]*this->components[k][j]*(k==0?1.:-1.);
    }
    return result;
  }

  LM<Field> transpose() {
    LM<Field> output;
    for(int mu1 = 0; mu1 < 4; mu1++) for(int mu2 = 0; mu2 < 4; mu2++) {
      output.components[mu1][mu2] = this->components[mu2][mu1];
    }
    return output;
  }

  Field trace() {
    Field output = 0;
    for(int i = 0; i < 4; i++) {
      output += (i==0?1.:-1.)*components[i][i];
    }
    return output;
  }
};


template <typename Field1, typename Field2>
LV<typename std::common_type<Field1, Field2>::type> operator/(const LV<Field1>& v, const Field2& scalar) {
  LV<typename std::common_type<Field1, Field2>::type> result;
  for (size_t i = 0; i < 4; ++i) {
    result.components[i] = v.components[i]/scalar;
  }
  return result;
}

template <typename Field1, typename Field2>
LV<typename std::common_type<Field1, Field2>::type> operator*(const Field1& scalar, const LV<Field2>& v) {
  LV<typename std::common_type<Field1, Field2>::type> output = v*scalar;
  return output;
}

template <class Field1, class Field2>
LV<typename std::common_type<Field1, Field2>::type> operator+(const LV<Field1>& v1, const LV<Field2>& v2) {
  LV<typename std::common_type<Field1, Field2>::type> output;
  for(int i = 0; i < 4; i++) {
    output.components[i] = v1.components[i] + v2.components[i];
  }
  return output;
}

template <class Field1, class Field2>
LV<typename std::common_type<Field1, Field2>::type> operator-(const LV<Field1>& v1, const LV<Field2>& v2) {
  LV<typename std::common_type<Field1, Field2>::type> output;
  for(int i = 0; i < 4; i++) {
    output.components[i] = v1.components[i] - v2.components[i];
  }
  return output;
}

template <class Field1, class Field2>
typename std::common_type<Field1, Field2>::type operator*(const LV<Field1>& v1, const LV<Field2>& v2) {
  typename std::common_type<Field1, Field2>::type output = 0;
  for(int i = 0; i < 4; i++) {
    output += (i==0?1.:-1.)*v1.components[i]*v2.components[i];
  }
  return output;
}

template <class Field1, class Field2>
LV<typename std::common_type<Field1, Field2>::type> operator*(const LM<Field1>& m, const LV<Field2>& v) {
  LV<typename std::common_type<Field1, Field2>::type> output = 0;
  for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) {
    output.components[i] += m.components[i][j]*v.components[j]*(j==0?1.:-1.);
  }
  return output;
}

#endif