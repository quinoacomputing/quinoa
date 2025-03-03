void ModUEqSolid::solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight,
double& dtMax, std::vector<double> &boundData) const
{
  Phase* vecPhase;
  double sL, sR;
  double rhoStar(0.), EStar(0.), stressTensor11Star(0.), stressTensor12Star(0.), stressTensor13Star(0.), vStar(0.), wStar(0.), pStar(0.);
  double uL = cellLeft.getMixture()->getVelocity().getX(), cL = 0.;
  double vL = cellLeft.getMixture()->getVelocity().getY(); double wL = cellLeft.getMixture()->getVelocity().getZ();
  double pL = cellLeft.getMixture()->getPressure(), rhoL = cellLeft.getMixture()->getDensity();
  double stressTensor11L = cellLeft.getMixture()->getStressTensor().getXX();
  double stressTensor12L = cellLeft.getMixture()->getStressTensor().getXY();
  double stressTensor13L = cellLeft.getMixture()->getStressTensor().getXZ();
  double uR = cellRight.getMixture()->getVelocity().getX(), cR = 0.;
  double vR = cellRight.getMixture()->getVelocity().getY(); double wR = cellRight.getMixture()->getVelocity().getZ();
  double pR = cellRight.getMixture()->getPressure(), rhoR = cellRight.getMixture()->getDensity();
  double stressTensor11R = cellRight.getMixture()->getStressTensor().getXX();
  double stressTensor12R = cellRight.getMixture()->getStressTensor().getXY();
  double stressTensor13R = cellRight.getMixture()->getStressTensor().getXZ();

  for (int k = 0; k < numberPhases; k++) {
  cL += cellLeft.getPhase(k)->getY() * cellLeft.getPhase(k)->getSquareLongitudinalWaveSpeed();
  cR += cellRight.getPhase(k)->getY() * cellRight.getPhase(k)->getSquareLongitudinalWaveSpeed();
  }

  cL = std::sqrt(cL);
  cR = std::sqrt(cR);
  cL *= 1.1; //To take into account the shear effect which is computationally heavy
  cR *= 1.1;

  //Low-Mach preconditioning
  double machRefMin(1.); //Default value without low-Mach preco.
  if (m_lowMach){
  lowMachSoundSpeed(machRefMin, uL, cL, uR, cR);
  }

  //Davies
  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);
  //For low-Mach (for general purpose machRefMin set to 1)
  if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
  if (std::fabs(sR)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxRight / std::fabs(sR));

  //Compute left and right mass flow rates and sM
  double mL = rhoL*(sL - uL), mR = rhoR*(sR - uR);
  double sM = (stressTensor11L - stressTensor11R + mL*uL - mR*uR) / (mL - mR);
  if (std::fabs(sM)<1.e-8) sM = 0.;

  //Compute the stress tensor components in the star region
  stressTensor11Star = (mR * stressTensor11L - mL * stressTensor11R + mL * mR * (uL - uR)) / (mR - mL);
  stressTensor12Star = (mR * stressTensor12L - mL * stressTensor12R + mL * mR * (vL - vR)) / (mR - mL);
  stressTensor13Star = (mR * stressTensor13L - mL * stressTensor13R + mL * mR * (wL - wR)) / (mR - mL);

  //Solution sampling
  double alpha, density, energy, pressure, lambda, lambdakStar, energyCompaction;
  if (sL >= 0.) {
    for (int k = 0; k < numberPhases; k++) {
    vecPhase = cellLeft.getPhase(k);
    alpha = vecPhase->getAlpha();
    density = vecPhase->getDensity();
    energy = vecPhase->getEnergy();
    lambda = vecPhase->getLambda();
    energyCompaction = vecPhase->getEnergyCompaction();
    static_cast<FluxUEqSolid*> (fluxBuff)->m_alpha[k] = alpha * sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_mass[k] = alpha * density * uL;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_energ[k] = alpha * density * (energy + energyCompaction) * uL;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_lambda[k] = lambda * uL;
    }
    for (int k = 0; k < numberSolids; k++) {
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setTensor(cellLeft.getPhase(k)->getCobase() * uL); //KS//S//Do a supersonic test case to know if uL,R or sM
    }
    double totalEnergy = cellLeft.getMixture()->getEnergy() + cellLeft.getMixture()->getEnergyCompaction() + cellLeft.getMixture()->getEnergyElastic() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setX(rhoL * uL * uL - stressTensor11L);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setY(rhoL * vL * uL - stressTensor12L);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setZ(rhoL * wL * uL - stressTensor13L);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_energMixture = (rhoL * totalEnergy - stressTensor11L) * uL - stressTensor12L * vL - stressTensor13L * wL;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_uStar = uL;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_vStar = vL;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_wStar = wL;
    // Boundary data for output
    boundData[VarBoundary::p] = pL;
    boundData[VarBoundary::rho] = rhoL;
    boundData[VarBoundary::velU] = uL;
    boundData[VarBoundary::velV] = vL;
    boundData[VarBoundary::velW] = wL;
  }

  else if (sR <= 0.) {
    for (int k = 0; k < numberPhases; k++) {
    vecPhase = cellRight.getPhase(k);
    alpha = vecPhase->getAlpha();
    density = vecPhase->getDensity();
    energy = vecPhase->getEnergy();
    lambda = vecPhase->getLambda();
    energyCompaction = vecPhase->getEnergyCompaction();
    static_cast<FluxUEqSolid*> (fluxBuff)->m_alpha[k] = alpha * sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_mass[k] = alpha * density * uR;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_energ[k] = alpha * density * (energy + energyCompaction) * uR;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_lambda[k] = lambda * uR;
    }
    for (int k = 0; k < numberSolids; k++) {
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setTensor(cellRight.getPhase(k)->getCobase() * uR); //KS//S//Do a supersonic test case to know if uL,R or sM
    }
    double totalEnergy = cellRight.getMixture()->getEnergy() + cellRight.getMixture()->getEnergyCompaction() + cellRight.getMixture()->getEnergyElastic() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setX(rhoR * uR * uR - stressTensor11R);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setY(rhoR * vR * uR - stressTensor12R);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setZ(rhoR * wR * uR - stressTensor13R);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_energMixture = (rhoR * totalEnergy - stressTensor11R) * uR - stressTensor12R * vR - stressTensor13R * wR;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_uStar = uR;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_vStar = vR;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_wStar = wR;
    // Boundary data for output
    boundData[VarBoundary::p] = pR;
    boundData[VarBoundary::rho] = rhoR;
    boundData[VarBoundary::velU] = uR;
    boundData[VarBoundary::velV] = vR;
    boundData[VarBoundary::velW] = wR;
  }

  else if (sM >= 0.) {
    //Compute left solution state
    double mkL;
    for (int k = 0; k < numberPhases; k++) {
    vecPhase = cellLeft.getPhase(k);
    alpha = vecPhase->getAlpha();
    density = vecPhase->getDensity();
    pressure = vecPhase->getPressure();
    mkL = density * (sL - uL);
    TB->rhokStar[k] = mkL / (sL - sM);
    TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
    TB->pkStar[k] = TB->eos[k]->computePressureIsentropic(pressure, density, TB->rhokStar[k]);
    // TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
    if (!m_lowMach) TB->ekStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]);
    else { TB->ekStar[k] = vecPhase->getEnergy(); }
    lambdakStar = vecPhase->getLambda() * (sL - uL) / (sL - sM);
    energyCompaction = vecPhase->getEnergyCompaction();
    static_cast<FluxUEqSolid*> (fluxBuff)->m_alpha[k] = alpha * sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_mass[k] = alpha * TB->rhokStar[k] * sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_energ[k] = alpha * TB->rhokStar[k] * (TB->ekStar[k] + energyCompaction) * sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_lambda[k] = lambdakStar * sM;
    }
    vStar = vL + (stressTensor12L - stressTensor12Star) / mL;
    wStar = wL + (stressTensor13L - stressTensor13Star) / mL;
    for (int k = 0; k < numberSolids; k++) {
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setTensor(cellLeft.getPhase(k)->getCobase() * sM);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setXX( (cellLeft.getPhase(k)->getCobase().getXX() * (uL - sL)
    + cellLeft.getPhase(k)->getCobase().getYX() * (vL - vStar)
    + cellLeft.getPhase(k)->getCobase().getZX() * (wL - wStar)) / (sM - sL) * sM);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setXY( (cellLeft.getPhase(k)->getCobase().getXY() * (uL - sL)
    + cellLeft.getPhase(k)->getCobase().getYY() * (vL - vStar)
    + cellLeft.getPhase(k)->getCobase().getZY() * (wL - wStar)) / (sM - sL) * sM);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setXZ( (cellLeft.getPhase(k)->getCobase().getXZ() * (uL - sL)
    + cellLeft.getPhase(k)->getCobase().getYZ() * (vL - vStar)
    + cellLeft.getPhase(k)->getCobase().getZZ() * (wL - wStar)) / (sM - sL) * sM);
    }
    double totalEnergy = cellLeft.getMixture()->getEnergy() + cellLeft.getMixture()->getEnergyCompaction() + cellLeft.getMixture()->getEnergyElastic() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();
    rhoStar = mL / (sL - sM);
    EStar = (-mL*totalEnergy - stressTensor11L*uL - stressTensor12L*vL - stressTensor13L*wL + stressTensor11Star*sM + stressTensor12Star*vStar + stressTensor13Star*wStar) / (rhoStar * (sM - sL));
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setX(rhoStar * sM * sM - stressTensor11Star);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setY(rhoStar * vStar * sM - stressTensor12Star);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setZ(rhoStar * wStar * sM - stressTensor13Star);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_energMixture = (rhoStar * EStar - stressTensor11Star) * sM - stressTensor12Star * vStar - stressTensor13Star * wStar;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_uStar = sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_vStar = vStar;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_wStar = wStar;
    // Boundary data for output
    pStar = mL*(sM - uL) + pL;
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vL;
    boundData[VarBoundary::velW] = wL;
  }

  else {
    //Compute right solution state
    double mkR;
    for (int k = 0; k < numberPhases; k++) {
    vecPhase = cellRight.getPhase(k);
    alpha = vecPhase->getAlpha();
    density = vecPhase->getDensity();
    pressure = vecPhase->getPressure();
    mkR = density * (sR - uR);
    TB->rhokStar[k] = mkR / (sR - sM);
    TB->eos[k]->verifyAndCorrectDensityMax(TB->rhokStar[k]);
    TB->pkStar[k] = TB->eos[k]->computePressureIsentropic(pressure, density, TB->rhokStar[k]);
    // TB->pkStar[k] = TB->eos[k]->computePressureHugoniot(pressure, density, TB->rhokStar[k]);
    if (!m_lowMach) TB->ekStar[k] = TB->eos[k]->computeEnergy(TB->rhokStar[k], TB->pkStar[k]);
    else { TB->ekStar[k] = vecPhase->getEnergy(); }
    lambdakStar = vecPhase->getLambda() * (sR - uR) / (sR - sM);
    energyCompaction = vecPhase->getEnergyCompaction();
    static_cast<FluxUEqSolid*> (fluxBuff)->m_alpha[k] = alpha * sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_mass[k] = alpha * TB->rhokStar[k] * sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_energ[k] = alpha * TB->rhokStar[k] * (TB->ekStar[k] + energyCompaction) * sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_lambda[k] = lambdakStar * sM;
    }
    vStar = vR + (stressTensor12R - stressTensor12Star) / mR;
    wStar = wR + (stressTensor13R - stressTensor13Star) / mR;
    for (int k = 0; k < numberSolids; k++) {
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setTensor(cellRight.getPhase(k)->getCobase() * sM);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setXX( (cellRight.getPhase(k)->getCobase().getXX() * (uR - sR)
    + cellRight.getPhase(k)->getCobase().getYX() * (vR - vStar)
    + cellRight.getPhase(k)->getCobase().getZX() * (wR - wStar)) / (sM - sR) * sM);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setXY( (cellRight.getPhase(k)->getCobase().getXY() * (uR - sR)
    + cellRight.getPhase(k)->getCobase().getYY() * (vR - vStar)
    + cellRight.getPhase(k)->getCobase().getZY() * (wR - wStar)) / (sM - sR) * sM);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setXZ( (cellRight.getPhase(k)->getCobase().getXZ() * (uR - sR)
    + cellRight.getPhase(k)->getCobase().getYZ() * (vR - vStar)
    + cellRight.getPhase(k)->getCobase().getZZ() * (wR - wStar)) / (sM - sR) * sM);
    }
    double totalEnergy = cellRight.getMixture()->getEnergy() + cellRight.getMixture()->getEnergyCompaction() + cellRight.getMixture()->getEnergyElastic() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();
    rhoStar = mR / (sR - sM);
    EStar = (-mR*totalEnergy - stressTensor11R*uR - stressTensor12R*vR - stressTensor13R*wR + stressTensor11Star*sM + stressTensor12Star*vStar + stressTensor13Star*wStar) / (rhoStar * (sM - sR));
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setX(rhoStar * sM * sM - stressTensor11Star);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setY(rhoStar * vStar * sM - stressTensor12Star);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setZ(rhoStar * wStar * sM - stressTensor13Star);
    static_cast<FluxUEqSolid*> (fluxBuff)->m_energMixture = (rhoStar * EStar - stressTensor11Star) * sM - stressTensor12Star * vStar - stressTensor13Star * wStar;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_uStar = sM;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_vStar = vStar;
    static_cast<FluxUEqSolid*> (fluxBuff)->m_wStar = wStar;
    //Boundary data for output
    pStar = mR*(sM - uR) + pR;
    boundData[VarBoundary::p] = pStar;
    boundData[VarBoundary::rho] = rhoStar;
    boundData[VarBoundary::velU] = sM;
    boundData[VarBoundary::velV] = vR;
    boundData[VarBoundary::velW] = wR;
  }

  //Contact discontinuity velocity
  static_cast<FluxUEqSolid*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions ****************
//****************************************************************************
void ModUEqSolid::solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double> &boundData) const
{
//Corresponds to a non-slipping wall
double sL(0.);
double stressTensor11Star(0.), stressTensor12Star(0.), stressTensor13Star(0.);
double uL = cellLeft.getMixture()->getVelocity().getX(), cL = 0.;
double vL = cellLeft.getMixture()->getVelocity().getY();
double wL = cellLeft.getMixture()->getVelocity().getZ(), rhoL = cellLeft.getMixture()->getDensity();
double stressTensor11L = cellLeft.getMixture()->getStressTensor().getXX();
double stressTensor12L = cellLeft.getMixture()->getStressTensor().getXY();
double stressTensor13L = cellLeft.getMixture()->getStressTensor().getXZ();
for (int k = 0; k < numberPhases; k++) {
cL += cellLeft.getPhase(k)->getY() * cellLeft.getPhase(k)->getSquareLongitudinalWaveSpeed();
}
cL = std::sqrt(cL);
cL *= 1.1; //To take into account the shear effect which is computationally heavy
//Low-Mach preconditioning
double machRefMin(1.); //Default value without low-Mach preco.
if (m_lowMach){
lowMachSoundSpeed(machRefMin, uL, cL);
}
//Davies
sL = std::min(uL - cL, -uL - cL);
//For low-Mach (for general purpose machRefMin set to 1)
if (std::fabs(sL)>1.e-3) dtMax = std::min(dtMax, machRefMin * dxLeft / std::fabs(sL));
//Compute the stress tensor components in the star region
stressTensor11Star = stressTensor11L - rhoL * (uL - sL) * uL;
stressTensor12Star = stressTensor12L - rhoL * (uL - sL) * vL;
stressTensor13Star = stressTensor13L - rhoL * (uL - sL) * wL;
//Solution
for (int k = 0; k < numberPhases; k++)
{
static_cast<FluxUEqSolid*> (fluxBuff)->m_alpha[k] = 0.;
static_cast<FluxUEqSolid*> (fluxBuff)->m_mass[k] = 0.;
static_cast<FluxUEqSolid*> (fluxBuff)->m_energ[k] = 0.;
static_cast<FluxUEqSolid*> (fluxBuff)->m_lambda[k] = 0.;
}
for (int k = 0; k < numberSolids; k++) {
static_cast<FluxUEqSolid*> (fluxBuff)->m_cobase[k].setTensor(0.);
}
static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setX(-stressTensor11Star);
static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setY(-stressTensor12Star);
static_cast<FluxUEqSolid*> (fluxBuff)->m_momentum.setZ(-stressTensor13Star);
static_cast<FluxUEqSolid*> (fluxBuff)->m_energMixture = 0.;
//Contact discontinuity velocity
static_cast<FluxUEqSolid*> (fluxBuff)->m_sM = 0.;
static_cast<FluxUEqSolid*> (fluxBuff)->m_uStar = 0.;
static_cast<FluxUEqSolid*> (fluxBuff)->m_vStar = 0.;
static_cast<FluxUEqSolid*> (fluxBuff)->m_wStar = 0.;
//Boundary data for output
boundData[VarBoundary::p] = stressTensor11Star;
boundData[VarBoundary::rho] = 0.;
boundData[VarBoundary::velU] = 0.;
boundData[VarBoundary::velV] = 0.;
boundData[VarBoundary::velW] = 0.;
}
