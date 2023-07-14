mkdir catkin_ws  //  创建catkin_ws的文件
cd catkin_ws    //打开catkin_ws这个文件
mkdir src      //然后在catkin_ws这个文件里再创建src这个文件

//将Uwb_Location放置在src文件夹里
cd catkin_ws
catkin_make        //在catkin_ws这个工作空间开始进行编译


//运行
//打开新的命令行
roscore
//再打开一个命令行
cd catkin_ws

source devel/setup.bash
rosrun Uwb_Location UwbLocation