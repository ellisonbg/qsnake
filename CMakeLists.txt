project(qsnake)

cmake_minimum_required(VERSION 2.6)

install(DIRECTORY qsnake/ DESTINATION lib/python/qsnake FILES_MATCHING PATTERN "*.py")
