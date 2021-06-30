################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/classes/KMC_Db.cpp \
../src/classes/KMC_Defs.cpp \
../src/classes/KMC_Kmer.cpp 

OBJS += \
./src/classes/KMC_Db.o \
./src/classes/KMC_Defs.o \
./src/classes/KMC_Kmer.o 

CPP_DEPS += \
./src/classes/KMC_Db.d \
./src/classes/KMC_Defs.d \
./src/classes/KMC_Kmer.d 


# Each subdirectory must supply rules for building sources it contributes
src/classes/%.o: ../src/classes/%.cpp src/classes/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


