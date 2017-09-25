      subroutine fillTable(table)
      use external_types
      implicit real (a-h,m,o-z)

      TYPE(FileGroupInfo),DIMENSION(11) :: table

      table(1)%initLink=10
      table(1)%tableType="cooling"
      table(1)%WDtype="DA"
      table(1)%flag=1
      table(1)%Z=0.001
      table(1)%ncol=7
      table(1)%nrow=650
      allocate(table(1)%ntrk(table(1)%ncol))
      allocate(table(1)%mass(table(1)%ncol))
      allocate(table(1)%prevTime(table(1)%ncol))
      allocate(table(1)%coolingTime(table(1)%ncol,table(1)%nrow))
      allocate(table(1)%luminosity(table(1)%ncol,table(1)%nrow))
      allocate(table(1)%effTemp(table(1)%ncol,table(1)%nrow))
      allocate(table(1)%gravAcc(table(1)%ncol,table(1)%nrow))

      table(2)%initLink=20
      table(2)%tableType="cooling"
      table(2)%WDtype="DA"
      table(2)%flag=2
      table(2)%Z=0.01
      table(2)%ncol=10
      table(2)%nrow=650
      allocate(table(2)%ntrk(table(2)%ncol))
      allocate(table(2)%mass(table(2)%ncol))
      allocate(table(2)%prevTime(table(2)%ncol))
      allocate(table(2)%coolingTime(table(2)%ncol,table(2)%nrow))
      allocate(table(2)%luminosity(table(2)%ncol,table(2)%nrow))
      allocate(table(2)%effTemp(table(2)%ncol,table(2)%nrow))
      allocate(table(2)%gravAcc(table(2)%ncol,table(2)%nrow))

      table(3)%initLink=30
      table(3)%tableType="cooling"
      table(3)%WDtype="DA"
      table(3)%flag=3
      table(3)%Z=0.03
      table(3)%ncol=8
      table(3)%nrow=650
      allocate(table(3)%ntrk(table(3)%ncol))
      allocate(table(3)%mass(table(3)%ncol))
      allocate(table(3)%prevTime(table(3)%ncol))
      allocate(table(3)%coolingTime(table(3)%ncol,table(3)%nrow))
      allocate(table(3)%luminosity(table(3)%ncol,table(3)%nrow))
      allocate(table(3)%effTemp(table(3)%ncol,table(3)%nrow))
      allocate(table(3)%gravAcc(table(3)%ncol,table(3)%nrow))

      table(4)%initLink=40
      table(4)%tableType="cooling"
      table(4)%WDtype="DA"
      table(4)%flag=3
      table(4)%Z=0.06
      table(4)%ncol=8
      table(4)%nrow=650
      allocate(table(4)%ntrk(table(4)%ncol))
      allocate(table(4)%mass(table(4)%ncol))
      allocate(table(4)%prevTime(table(4)%ncol))
      allocate(table(4)%coolingTime(table(4)%ncol,table(4)%nrow))
      allocate(table(4)%luminosity(table(4)%ncol,table(4)%nrow))
      allocate(table(4)%effTemp(table(4)%ncol,table(4)%nrow))
      allocate(table(4)%gravAcc(table(4)%ncol,table(4)%nrow))

      table(5)%initLink=90
      table(5)%tableType="cooling"
      table(5)%WDtype="DB"
      table(5)%flag=1
      table(5)%Z=0.001
      table(5)%ncol=7
      table(5)%nrow=400
      allocate(table(5)%ntrk(table(5)%ncol))
      allocate(table(5)%mass(table(5)%ncol))
      allocate(table(5)%prevTime(table(5)%ncol))
      allocate(table(5)%coolingTime(table(5)%ncol,table(5)%nrow))
      allocate(table(5)%luminosity(table(5)%ncol,table(5)%nrow))
      allocate(table(5)%effTemp(table(5)%ncol,table(5)%nrow))
      allocate(table(5)%gravAcc(table(5)%ncol,table(5)%nrow))

      table(6)%initLink=100
      table(6)%tableType="cooling"
      table(6)%WDtype="DB"
      table(6)%flag=2
      table(6)%Z=0.01
      table(6)%ncol=9
      table(6)%nrow=400
      allocate(table(6)%ntrk(table(6)%ncol))
      allocate(table(6)%mass(table(6)%ncol))
      allocate(table(6)%prevTime(table(6)%ncol))
      allocate(table(6)%coolingTime(table(6)%ncol,table(6)%nrow))
      allocate(table(6)%luminosity(table(6)%ncol,table(6)%nrow))
      allocate(table(6)%effTemp(table(6)%ncol,table(6)%nrow))
      allocate(table(6)%gravAcc(table(6)%ncol,table(6)%nrow))

      table(7)%initLink=110
      table(7)%tableType="cooling"
      table(7)%WDtype="DB"
      table(7)%flag=3
      table(7)%Z=0.06
      table(7)%ncol=9
      table(7)%nrow=400
      allocate(table(7)%ntrk(table(7)%ncol))
      allocate(table(7)%mass(table(7)%ncol))
      allocate(table(7)%prevTime(table(7)%ncol))
      allocate(table(7)%coolingTime(table(7)%ncol,table(7)%nrow))
      allocate(table(7)%luminosity(table(7)%ncol,table(7)%nrow))
      allocate(table(7)%effTemp(table(7)%ncol,table(7)%nrow))
      allocate(table(7)%gravAcc(table(7)%ncol,table(7)%nrow))

      table(8)%initLink=120
      table(8)%tableType="color"
      table(8)%WDtype="ONe"

      table(9)%initLink=127
      table(9)%tableType="cooling"
      table(9)%WDtype="ONe"

      table(10)%initLink=131
      table(10)%tableType="colors"
      table(10)%WDtype="DB"
      table(10)%ncol=7
      table(10)%nrow=60
      allocate(table(10)%ntrk(table(10)%ncol))
      allocate(table(10)%mass(table(10)%ncol))
      allocate(table(10)%luminosity(table(10)%ncol,table(10)%nrow))
      allocate(table(10)%color_U(table(10)%ncol,table(10)%nrow))
      allocate(table(10)%color_B(table(10)%ncol,table(10)%nrow))
      allocate(table(10)%color_V(table(10)%ncol,table(10)%nrow))
      allocate(table(10)%color_R(table(10)%ncol,table(10)%nrow))
      allocate(table(10)%color_I(table(10)%ncol,table(10)%nrow))
      allocate(table(10)%color_J(table(10)%ncol,table(10)%nrow))

      table(11)%initLink=60
      table(11)%tableType="colors"
      table(11)%WDtype="DA"
      table(11)%ncol=10
      table(11)%nrow=650
      allocate(table(11)%ntrk(table(11)%ncol))
      allocate(table(11)%mass(table(11)%ncol))
      allocate(table(11)%luminosity(table(11)%ncol,table(11)%nrow))
      allocate(table(11)%color_U(table(11)%ncol,table(11)%nrow))
      allocate(table(11)%color_B(table(11)%ncol,table(11)%nrow))
      allocate(table(11)%color_V(table(11)%ncol,table(11)%nrow))
      allocate(table(11)%color_R(table(11)%ncol,table(11)%nrow))
      allocate(table(11)%color_I(table(11)%ncol,table(11)%nrow))
      allocate(table(11)%color_J(table(11)%ncol,table(11)%nrow))

      return
      end